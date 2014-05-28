all:

html:
	ipython nbconvert --to html *.ipynb

markdown:
	ipython nbconvert --to markdown *.ipynb