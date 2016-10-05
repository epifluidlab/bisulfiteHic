#!/bin/bash
cd jbwa/jbwa-1.0.0/ && make && cd ../../ && mvn compile assembly:single

if [ -z $BISHIC ]
then
	export BISHIC=${PWD}
	echo "export BISHIC=${BISHIC}" >> ~/.bash_profile
	source ~/.bash_profile
fi

chmod 755 bisulfitehicMap
cd ./target && ln -s $(ls -t | head -n1) bisulfitehic-default.jar
