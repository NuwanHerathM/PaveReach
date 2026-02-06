all: build_docker chmod

build_docker:
	docker build -t globalqe .

chmod:
	chmod +x Fig*.sh Tab*.sh

run_Fig5:
	docker run -v .:/app globalqe examples/ex_5-2_stability_controller.jl 0.1 --with_luxor

clean:
	rm -f *.png *.log