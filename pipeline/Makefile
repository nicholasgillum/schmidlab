all: sweave
	cd report; make $(MKFLGS)

sweave: force_look
	mv preprocessing.tex report/preprocessing.tex;
	
force_look:
	true

clean:
	cd tmp; rm -rf preprocessing.*
	rm -f preprocessing.tex
	rm -f tmp/*
