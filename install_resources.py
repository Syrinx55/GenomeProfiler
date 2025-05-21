from configparser import SectionProxy
from pycurl import Curl
import requests
import tarfile
from tempfile import TemporaryFile
from lxml import html


URL_PLSDB_META = "https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_meta.tar.gz"
URL_MOBILEOG_DB = "https://mobileogdb.flsi.cloud.vt.edu/entries/database_download"


# TODO test
def install_plsdb_meta(to_dir: str):
    response = requests.get(URL_PLSDB_META, stream=True)
    if response.status_code != 200:
        raise ConnectionError(
            f"Request to {URL_PLSDB_META} returned code {response.status_code}"
        )

    with TemporaryFile() as f_download:
        f_download.write(response.raw.read())

        with tarfile.open(fileobj=f_download) as f_tar:
            f_tar.extractall(path=to_dir)


def install_mobileog_db(to_dir: str):
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

    with TemporaryFile() as f_download:
        cl = Curl()
        cl.setopt(cl.URL, targz_url)
        cl.setopt(cl.WRITEDATA, f_download)
        cl.perform()
        cl.close()

        with tarfile.open(fileobj=f_download) as f_tar:
            f_tar.extractall(path=to_dir)


def install_tncentral_db(to_dir: str):
    pass


def install_resources(config: SectionProxy):
    print("Installing PLSDB meta...")
    install_plsdb_meta(config["plsdb_meta_dir"])
    print("Installed PLSDB meta.")

    print("Installing mobileOG-db...")
    # FIXME reference instead of hard-code path
    install_mobileog_db("mobileOG-db-beatrix-1.6")
    print("Installed mobileOG-db.")

    print("Installing TnCentral database...")
    install_tncentral_db(config["tncentral_db"])
    print("Installed TnCentral database.")
