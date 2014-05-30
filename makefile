all: html markdown

html:
	ipython nbconvert --to html *.ipynb

markdown:
	ipython nbconvert --to markdown *.ipynb

pdf:
	ipython nbconvert --to latext --post PDF *.ipynb
