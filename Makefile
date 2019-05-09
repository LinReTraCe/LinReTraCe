all: binR linretrace make_include

binR:
	if [ ! -d bin ] ; then mkdir bin ; fi

linretrace: make_include
	cd src_linretrace/; make

clean:
	cd src_linretrace/; make clean
	cd src_linretrace/digamma; make clean

pristine:
	cd src_linretrace/; make clean
	cd src_linretrace/digamma; make clean
	rmdir bin

ctags:
	ctags -R --exclude=src_linretrace/digamma --exclude=src_pp --exclude=src_linretrace/deprecated

