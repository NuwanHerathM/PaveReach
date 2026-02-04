build:
	docker build -t globalqe .

run:
	docker run -v .:/app globalqe --with_luxor

clean:
	rm -f *.png *.log