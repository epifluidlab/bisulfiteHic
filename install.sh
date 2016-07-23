#!/bin/bash
cd jbwa/jbwa-1.0.0/ && make && cd ../../ && sudo mvn compile assembly:single

if [ -z $BISHIC ]
then
	export BISHIC=${PWD}
	echo "export BISHIC=${BISHIC}" >> ~/.bash_profile
	source ~/.bash_profile
fi

chmod 755 bisulfitehicMap