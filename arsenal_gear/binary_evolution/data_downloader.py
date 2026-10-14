"""
data_downloader.py
================================

This file defines the interface to download output from various binary
evolution models.
"""

from pathlib import Path

from arsenal_gear.utils import downloader, untar, unzip


class BPASSDownloader:
    """
    Downloads and extracts the BPASS stellar models needed for conversion.
    """

    # basic options for BPASS stellar models
    # bpass_url = "https://www.dropbox.com/scl/fo/mpuas1xh5owmdadu0vpev/h?dl=0" + \
    #            "&e=1&rlkey=7vlk7ra6kvoztzmae8wr34kmz"
    # Change dl=0 to dl=1 to force download
    bpass_url = "https://www.dropbox.com/scl/fo/mpuas1xh5owmdadu0vpev/h?dl=1&e=1&rlkey=7vlk7ra6kvoztzmae8wr34kmz"
    # name of the release on the server, and the only member of it that we need
    zip_name = "BPASSv2.2.1_release-07-18-Tuatara.zip"
    tar_name = "bpass-v2.2-newmodels.tar.gz"
    # directory created by extracting tar_name
    models_name = "NEWBINMODS"

    def __init__(
        self,
        bpass_dir: str | None = None,
        force_download: bool = False,
        delete_zip: bool = True,
    ) -> None:
        """
        Args:
            bpass_dir: the directory for the BPASS models
                       (default: binary_evolution/data/BPASS)
            force_download: download the release if no local copy is found
            delete_zip: delete the (~13 GB) release zip once the models
                        tar file has been extracted from it
        """

        self.force_download = force_download
        self.delete_zip = delete_zip
        # make data storage directory if it doesn't exist
        if bpass_dir is None:
            bpass_dir = Path(__file__).parent.absolute() / "data/BPASS"
            bpass_dir.mkdir(parents=True, exist_ok=True)
        self.dir = Path(bpass_dir)

    def download(self, message="Downloading stellar models...") -> None:
        """
        Method for downloading the BPASS release from the web.

        Args:
            message (str): Optional message to display before downloading.

        Raises:
            Exception: If the download fails.

        """

        self.dir.mkdir(parents=True, exist_ok=True)
        downloader(self.dir / self.zip_name, self.bpass_url, message)

    def get_stellar_models(self) -> None:
        """
        Look for the stellar models for the stellar population and
        download/unzip/untar them if they are not available.
        """
        models_dir = self.dir / self.models_name
        tar_path = self.dir / self.tar_name
        zip_path = self.dir / self.zip_name

        if models_dir.is_dir():
            print("BPASS data already available. Ready to start converting.")
            return
        if not tar_path.is_file():
            if not zip_path.is_file():
                if not self.force_download:
                    print(
                        f"BPASS data not found in {self.dir}, "
                        "set force_download = True to download"
                    )
                    return
                self.download()
            unzip(
                zip_path, self.dir, target_file=self.tar_name, delete_zip=self.delete_zip
            )
        untar(tar_path, self.dir, delete_tar=self.delete_zip)
        if not models_dir.is_dir():
            raise FileNotFoundError(
                f"Extracting {tar_path} did not produce {models_dir}. "
                "The layout of the BPASS release may have changed."
            )
        print("BPASS data now available. Ready to start converting.")


class MPADownloader:
    """
    Reads in MPA/Bonn stellar model files for use with a discrete stellar population.
    """

    mpa_url = "https://wwwmpa.mpa-garching.mpg.de/stellgrid/data/BonnGrids/"

    def __init__(
        self,
        mpa_dir: str | None = None,
        username: str | None = None,
        password: str | None = None,
        force_download: bool = False,
    ) -> None:
        """
        Args:
            mpa_dir: the directory for the MPA/Bonn models
                     (default: binary_evolution/data/MPA)
            username, password: credentials for the MPA/Bonn data server,
                                only needed to download

        """

        self.force_download = force_download
        self.username = username
        self.password = password
        # make data storage directory if it doesn't exist
        if mpa_dir is None:
            mpa_dir = Path(__file__).parent.absolute() / "data/MPA"
            mpa_dir.mkdir(parents=True, exist_ok=True)
        self.dir = Path(mpa_dir)

    def download(self, url=mpa_url, message="Downloading stellar models...") -> None:
        """
        Method for downloading MPA/Bonn stellar models from the web.

        Args:
            url (str): The URL of the MPA/Bonn repository.
            message (str): Optional message to display before downloading.

        Raises:
            Exception: If the download fails.

        """

        model_sets = [
            "MW/primary",
            "MW/secondary",
            "single_MW/single_MW",
            "LMC/primary",
            "LMC/secondary",
            "single_LMC/single_LMC",
            "SMC/primary",
            "SMC/secondary",
            "single_SMC/single_SMC",
        ]

        for model_set in model_sets:

            download_url = url + model_set + ".tar"
            fname = self.dir / f"{model_set}.tar"
            fname.parent.mkdir(parents=True, exist_ok=True)

            downloader(fname, download_url, message, self.username, self.password)
