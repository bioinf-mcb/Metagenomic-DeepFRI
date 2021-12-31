import pathlib
import ftplib
import multiprocessing
import requests
import shutil


url = "ftp.wwpdb.org"
ftp = ftplib.FTP(url)
ftp.login("anon", "anon")

ftp.cwd("/pub/pdb/data/structures/all/mmCIF/")
files_to_download = ftp.nlst()
files_to_download = files_to_download[2:]

download_path = pathlib.Path("./PDB_Compressed")
download_path.mkdir(exist_ok=True)
files_downloaded = [x.name for x in list(download_path.glob("*cif.gz"))]

files_to_download = [x for x in (set(files_to_download) - set(files_downloaded))]


def download_file(name):
    with requests.get("https://files.rcsb.org/download/"+name, stream=True) as r:
        with open(download_path / name, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


with multiprocessing.Pool(processes=64) as p:
    p.map(download_file, files_to_download)
