"""
scraper.py
==========

This file contains various utility functions for scraping data from online sources
Functions:
    downloader: Method for downloading data from the web
"""

import requests
from tqdm import tqdm

def downloader(fname, url, message):
    """
    Method for downloading isochrone data from the web where available.

    Args:
        fname (str or Path): The file name or path to save the downloaded file.
        url (str): The URL of the file to download.
        message (str): Optional message to display before downloading.

    Raises:
        Exception: If the download fails.

    """
    if message is not None:
        print(message)

    try:
        response = requests.get(url, stream=True, timeout=10)
    except requests.exceptions.Timeout as e:
        raise TimeoutError('Request timed out. Check internet connection.') from e
    except requests.exceptions.ConnectionError as e:
        raise ConnectionError('Connection error. Check internet connection.') from e
    except requests.exceptions.HTTPError as e:
        raise RuntimeError(f'HTTP error occurred: {e}') from e
    except requests.exceptions.TooManyRedirects as e:
        raise RuntimeError('Too many redirects. Check the URL.') from e
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f'Download failed: {e}') from e

    # Get file size
    total_size = int(response.headers.get('content-length', 0))
    # create a progress bar
    tqdm_args = {"desc": "Downloading", "total": total_size, "unit": 'B',
                    "unit_scale": True, "unit_divisor": 1024}
    # write the file
    with open(fname, 'wb') as f, tqdm(**tqdm_args) as prog_bar:
        for chunk in response.iter_content(chunk_size=1024):
            f.write(chunk)
            prog_bar.update(len(chunk))
