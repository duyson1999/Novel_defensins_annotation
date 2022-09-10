import os


def read_txt(path):
	return open(path, 'r').read()


def write_txt(path, txt):
	return open(path, 'w').write(txt)


def mkdir(path):
	return os.makedirs(path, exist_ok=True)


def check_size(path):
	return os.stat(path).st_size