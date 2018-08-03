all: binR linretrace make_include

binR:
	if [ ! -d bin ] ; then mkdir bin ; fi

linretrace: make_include
	cd src/; make

clean:
	cd src/; make clean

pristine:
	cd src/; make clean
	rmdir bin

ctags:
	ctags -R --exclude=src/digamma
