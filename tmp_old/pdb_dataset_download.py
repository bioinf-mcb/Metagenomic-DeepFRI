import pathlib
import ftplib
import multiprocessing
import requests


def chunkify(lst, n):
    return [lst[i::n] for i in range(n)]


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

chunks = chunkify(files_to_download, 1024)


def download_files(file_list):
    # ftp download speed is unfeasible
    for name in file_list:
        r = requests.get("https://files.rcsb.org/download/"+name, allow_redirects=True)
        with open(download_path / name, 'wb') as f:
            f.write(r.content)


with multiprocessing.Pool(processes=64) as p:
    p.map(download_files, chunks)
