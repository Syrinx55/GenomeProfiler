from configparser import SectionProxy
import requests
import tarfile
from tempfile import TemporaryFile


URL_PLSDB_META = "https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_meta.tar.gz"


# TODO test
def install_plsdb_meta(to_dir: str):
    with TemporaryFile() as f_download:
        response = requests.get(URL_PLSDB_META, stream=True)
        if response.status_code != 200:
            raise ConnectionError(
                f"Request to {URL_PLSDB_META} returned code {response.status_code}"
            )

        f_download.write(response.raw.read())

        with tarfile.open(fileobj=f_download) as f_tar:
            f_tar.extractall(path=to_dir)


def install_mobileog_db(to_dir: str):
    pass


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
