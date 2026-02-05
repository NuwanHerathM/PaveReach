build:
	docker build -t globalqe .

run_Fig5:
	docker run -v .:/app globalqe examples/ex_5-2_stability_controller.jl 0.1 --with_luxor


clean:
	rm -f *.png *.log