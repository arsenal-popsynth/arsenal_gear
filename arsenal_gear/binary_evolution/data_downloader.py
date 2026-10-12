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
    Reads in BPASS stellar model files for use with a discrete stellar population.
    """

    # basic options for BPASS stellar models
    # bpass_url = "https://www.dropbox.com/scl/fo/mpuas1xh5owmdadu0vpev/h?dl=0" + \
    #            "&e=1&rlkey=7vlk7ra6kvoztzmae8wr34kmz"
    # Change dl=0 to dl=1 to force download
    bpass_url = "https://www.dropbox.com/scl/fo/mpuas1xh5owmdadu0vpev/h?dl=1&e=1&rlkey=7vlk7ra6kvoztzmae8wr34kmz"

    def __init__(
        self,
        bpass_dir: str,
        force_download: bool = False,
    ) -> None:
        """
        Args:
            bpass_dir: the directory for the BPASS models

        """

        self.force_download = force_download
        self.dir = Path(bpass_dir)

        super().__init__()

    def download(self, url=bpass_url, message="Downloading stellar models...") -> None:
        """
        Method for downloading BPASS stellar models from the web.

        Args:
            url (str): The URL of the BPASS repository.
            message (str): Optional message to display before downloading.

        Raises:
            Exception: If the download fails.

        """

        self.dir.mkdir(parents=True, exist_ok=True)
        fname = self.dir / "bpass_v2.2.zip"

        downloader(fname, url, message)

    def get_stellar_models(
        self,
        tar_name="bpass-v2.2-newmodels.tar.gz",
        zip_name="bpass_v2.2.zip",
        url=bpass_url,
    ) -> tuple:
        """
        Look for the stellar models for the stellar population
        and download/unzip/untar them if they are not available.
        """
        if (self.dir / "NEWBINMODS").is_dir():
            print("BPASS data already available. Ready to start converting.")
        else:
            print("BPASS data not available at", self.dir / "NEWBINMODS")
            print("Looking for a tar file...")
            if (self.dir / tar_name).is_file():
                untar(self.dir / tar_name, self.dir, delete_tar=False)
                print("BPASS data now available. Ready to start converting.")

            else:
                print("tar file not available at", self.dir)
                print("Looking for a zip file...")
                if (self.dir / zip_name).is_file():
                    unzip(
                        self.dir / zip_name,
                        self.dir,
                        target_file=tar_name,
                        delete_zip=False,
                        inspect=False,
                    )
                    untar(self.dir / tar_name, self.dir, delete_tar=False)
                    print("BPASS data now available. Ready to start converting.")

                else:
                    print("zip file not available at", self.dir)
                    if self.force_download:
                        self.download(url)
                        unzip(
                            self.dir / zip_name,
                            self.dir,
                            target_file=tar_name,
                            delete_zip=False,
                            inspect=False,
                        )
                        untar(self.dir / tar_name, self.dir, delete_tar=False)
                        print("BPASS data now available. Ready to start converting.")

                    else:
                        print("Set force_download = True to download the files")

        return


class MPADownloader:
    """
    Reads in MPA/Bonn stellar model files for use with a discrete stellar population.
    """

    mpa_url = "https://wwwmpa.mpa-garching.mpg.de/stellgrid/data/BonnGrids/"

    def __init__(
        self,
        mpa_dir: str,
        username: str,
        password: str,
        force_download: bool = False,
    ) -> None:
        """
        Args:
            mpa_dir: the directory for the MPA/Bonn models

        """

        self.force_download = force_download
        self.username = username
        self.password = password
        self.dir = Path(mpa_dir)

        super().__init__()

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
