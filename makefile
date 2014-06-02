all: html markdown

tar:
	tar -zcvf data_1.tar.gz data_1
	tar -zcvf data_2.tar.gz data_2
	tar -zcvf data_3.tar.gz data_3
	tar -zcvf data_4.tar.gz data_4

untar:
	tar -zxvf data_1.tar.gz
	tar -zxvf data_2.tar.gz
	tar -zxvf data_3.tar.gz
	tar -zxvf data_4.tar.gz
	mkdir data
	mv data_1/* data
	mv data_2/* data
	mv data_3/* data
	mv data_4/* data

html:
	ipython nbconvert --to html *.ipynb

markdown:
	ipython nbconvert --to markdown *.ipynb

pdf:
	ipython nbconvert --to latext --post PDF *.ipynb
