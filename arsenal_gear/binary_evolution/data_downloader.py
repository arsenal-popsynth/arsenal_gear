"""
data_converter.py
================================

This file defines the interface to download output from various binary
evolution models.
"""

import os
import tarfile
from zipfile import ZipFile
from pathlib import Path
import requests
from tqdm import tqdm

class download_BPASS_models:
    """
    Reads in BPASS stellar model files for use with a discrete stellar population.
    """

    # basic options for BPASS stellar models
    # bpass_url = "https://www.dropbox.com/scl/fo/mpuas1xh5owmdadu0vpev/h?dl=0" + \
    #            "&e=1&rlkey=7vlk7ra6kvoztzmae8wr34kmz"
    # Change dl=0 to dl=1 to force download
    bpass_url = (
        "https://www.dropbox.com/scl/fo/mpuas1xh5owmdadu0vpev/h?dl=1"
        + "&e=1&rlkey=7vlk7ra6kvoztzmae8wr34kmz"
    )

    def __init__(
        self,
        bpass_dir: str,
        force_download: bool = False,
    ) -> None:
        """
        Args:
            bpass_dir: the directory for the BPASS models

        """

        self.download = force_download

        if bpass_dir[-1] == "/":
            self.dir: str = bpass_dir
        else:
            self.dir: str = bpass_dir + "/"

        super().__init__()

    def downloader(
        self, url=bpass_url, message="Downloading stellar models..."
    ) -> None:
        """
        Method for downloading BPASS stellar models from the web.

        Args:
            url (str): The URL of the BPASS repository.
            message (str): Optional message to display before downloading.

        Raises:
            Exception: If the download fails.

        """
        if message is not None:
            print(message)

        download_url = url
        if self.dir[-1] == "/":
            fname = self.dir + "bpass_v2.2.zip"
        else:
            fname = self.dir + "/bpass_v2.2.zip"

        try:
            response = requests.get(download_url, stream=True, timeout=10)
        except requests.exceptions.Timeout as e:
            raise requests.exceptions.Timeout(
                f"Request timed out. Check your internet connection. {e}"
            )
        except requests.exceptions.RequestException as e:
            raise requests.exceptions.RequestException(f"Request failed: {e}")

        # Get file size
        total_size = int(response.headers.get("content-length", 0))
        # create a progress bar
        tqdm_args = {'desc': 'Downloading', 'total': total_size, 'unit': 'B',
                        'unit_scale': True, 'unit_divisor': 1024}
        # write the file
        with open(fname, "wb") as f, tqdm(**tqdm_args) as prog_bar:
            for chunk in response.iter_content(chunk_size=1024):
                f.write(chunk)
                prog_bar.update(len(chunk))

    def unzip(
        self, zip_name: str, target_file=None, delete_zip=False, inspect=False
    ) -> None:
        """
        Un-compress the zip file downloaded from the BPASS dropbox.

        Args:
            zip_name    (str): Name of expected zip file
            target_file (str): Name of expected tar file in zip file
                               If none, extract all tar files
            delete_zip (bool): Delete the zip file after extracting the tar file
            inspect    (bool): Print the names of the tar files in the zip file
        """
        if self.dir[-1] == "/":
            fname = self.dir + zip_name
        else:
            fname = self.dir + "/" + zip_name

        with ZipFile(fname, "r") as zip_archive:
            file_names = zip_archive.namelist()
            if inspect:
                for name in file_names:
                    print(name)
            else:
                if target_file in file_names:

                    zip_archive.extract(target_file, path=self.dir)
                elif target_file is None:
                    zip_archive.extractall(path=self.dir)

        if delete_zip and fname.exists():
            print("Deleting", fname)
            fname.unlink()

    def untar(self, tar_name: str, delete_tar=False) -> None:
        """
        Un-compress the tar file extracted from the BPASS zip file.

        Args:
            tar_name    (str): Name of expected tar file
            delete_tar (bool): Delete the tar file after extracting
        """
        if self.dir[-1] == "/":
            fname = self.dir + tar_name
        else:
            fname = self.dir + "/" + tar_name

        if not tarfile.is_tarfile(fname):
            raise OSError(
                f"{fname} is not a valid tar.gz file. "
                "Try again with `force_download=True`"
            )
        with tarfile.open(fname, "r:gz") as tar:
            print(f"Extracting {fname}...")
            tar.extractall(path=self.dir)

        if delete_tar and fname.exists():
            print("Deleting", fname)
            fname.unlink()


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
        if os.path.isdir(self.dir + "NEWBINMODS"):
            print("BPASS data already available. Ready to start converting.")
        else:
            print("BPASS data not available at", self.dir + "NEWBINMODS")
            print("Looking for a tar file...")
            if os.path.isfile(self.dir + tar_name):
                self.untar(tar_name, delete_tar=False)
                print("BPASS data now available. Ready to start converting.")

            else:
                print("tar file not available at", self.dir)
                print("Looking for a zip file...")
                if os.path.isfile(self.dir + zip_name):
                    self.unzip(
                        zip_name, target_file=tar_name, delete_zip=False, inspect=False
                    )
                    self.untar(tar_name, delete_tar=False)
                    print("BPASS data now available. Ready to start converting.")

                else:
                    print("zip file not available at", self.dir)
                    if self.download:
                        self.downloader(url)
                        self.unzip(
                            zip_name,
                            target_file=tar_name,
                            delete_zip=False,
                            inspect=False,
                        )
                        self.untar(tar_name, delete_tar=False)
                        print("BPASS data now available. Ready to start converting.")

                    else:
                        print("Set force_download = True to download the files")

        return 

        