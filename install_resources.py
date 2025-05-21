from configparser import SectionProxy
from lxml import html
from pathlib import Path
from pycurl import Curl
import requests
import tarfile
from tempfile import TemporaryFile
from zipfile import ZipFile


URL_PLSDB_META = "https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_meta.tar.gz"
URL_MOBILEOG_DB = "https://mobileogdb.flsi.cloud.vt.edu/entries/database_download"
URL_TNCENTRAL_DB = "https://tncentral.ncc.unesp.br/api/download_blast/nc/tn_in_is"

# Keys of config whose values should be files that exist
PLSDB_META_VSBF = [
    "plsdb_sketch_path",
    "assembly_csv",
    "nuccore_csv",
    "plasmidfinder.csv",
    "typing.csv",
    "taxonomy.csv",
    "biosample.csv",
    "changes.tsv",
    "amr.tsv",
    "plsdb_mashdb_sim.tsv",
]

MOBILEOG_DB_VSBF = [
    "mobileog_db_faa",
    "mobileog_db_csv",
]

TNCENTRAL_DB_VSBF = [
    "tncentral_fasta",
]


def _files_not_installed(config: SectionProxy, values_should_be_files: list[str]) -> list[Path]:
    not_installed = []

    for key in values_should_be_files:
        path = Path(config[key])
        if not path.is_file():
            not_installed.append(path)

    return not_installed


def _curl_download(fileobj, url: str):
    cl = Curl()
    cl.setopt(cl.URL, url)
    cl.setopt(cl.WRITEDATA, fileobj)
    cl.perform()
    cl.close()


def _download_targz_and_extract(to_dir: str, url: str):
    with TemporaryFile() as f_download:
        _curl_download(f_download, url)

        with tarfile.open(fileobj=f_download) as f_tar:
            f_tar.extractall(path=to_dir)


def _download_zip_and_extract(to_dir: str, url: str):
    with TemporaryFile() as f_download:
        _curl_download(f_download, url)

        with ZipFile(f_download) as f_zip:
            f_zip.extractall(path=to_dir)


def _install_plsdb_meta(to_dir: str):
    _download_targz_and_extract(to_dir, URL_PLSDB_META)

def _install_mobileog_db(to_dir: str):
    targz_url = None
    with requests.get(URL_MOBILEOG_DB) as page_response:
        if page_response.status_code != 200:
            raise ConnectionError(
                f"Request to {URL_MOBILEOG_DB} returned code {page_response.status_code}"
            )
        tree_page = html.fromstring(page_response.content)
        targz_url = tree_page.xpath(
            "(/html/body/main/ul/li[1]/div/div[@class='btn-group'])[5]/a/@href"
        )[0]
    _download_zip_and_extract(to_dir, targz_url)


def _install_tncentral_db(to_dir: str):
    _download_zip_and_extract(to_dir, URL_TNCENTRAL_DB)


def install_plsdb_meta_if_not_present(to_dir: str, config: SectionProxy):
    not_installed = _files_not_installed(config, PLSDB_META_VSBF)
    if len(not_installed) == 0:
        print("(PLSDB meta already installed.)")
        return
    _install_plsdb_meta(to_dir)


def install_mobileog_db_if_not_present(to_dir: str, config: SectionProxy):
    not_installed = _files_not_installed(config, MOBILEOG_DB_VSBF)
    if len(not_installed) == 0:
        print("(mobileOG-db already installed.)")
        return
    _install_mobileog_db(to_dir)


def install_tncentral_db_if_not_present(to_dir: str, config: SectionProxy):
    not_installed = _files_not_installed(config, TNCENTRAL_DB_VSBF)
    if len(not_installed) == 0:
        print("(TnCentral database already installed.)")
        return
    _install_tncentral_db(to_dir)


def install_resources(config: SectionProxy, resource_path: str):
    print("Installing PLSDB meta...")
    install_plsdb_meta_if_not_present(config["plsdb_meta_dir"], config)
    print("Installed PLSDB meta.")

    print("Installing mobileOG-db...")
    install_mobileog_db_if_not_present(resource_path + "/mobileOG-db-beatrix-1.6", config)
    print("Installed mobileOG-db.")

    print("Installing TnCentral database...")
    install_tncentral_db_if_not_present(config["tncentral_db"], config)
    print("Installed TnCentral database.")
