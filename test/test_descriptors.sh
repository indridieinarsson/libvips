#!/bin/sh

# test the various restartable loaders

# set -x
set -e

. ./variables.sh

if test_supported jpegload; then
	./test_descriptors $image
fi

if test_supported heifload; then
	./test_descriptors $test_images/Example1.heic
fi

if test_supported gifload; then
	./test_descriptors $test_images/cogs.gif
fi

if test_supported pdfload; then
	./test_descriptors $test_images/ISO_12233-reschart.pdf 
fi

if test_supported pngload; then
	./test_descriptors $test_images/sample.png
fi

if test_supported webpload; then
	./test_descriptors $test_images/1.webp
fi
