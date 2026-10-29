all: build_docker decompress chmod

build_docker:
	docker build -t pavereach .

decompress:
	if [ ! -f examples/models/normed_etcs_2x250_nodes.onnx ]; then \
		gzip -d examples/models/normed_etcs_2x250_nodes.onnx.gz; \
	fi
	if [ ! -f examples/models/ACASXU_run2a_1_1_batch_2000.nnet ]; then \
		gzip -d examples/models/ACASXU_run2a_1_1_batch_2000.nnet.gz; \
	fi

chmod:
	chmod +x Figure*.sh Table1.sh

run_Figure5a:
	docker run -v .:/app pavereach examples/exp_5-1_global_robustness_landscape_class_1.jl

clean:
	rm -r examples/outputs