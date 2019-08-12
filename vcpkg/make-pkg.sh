#!/bin/sh

if [ -z "$1" ]
then
	echo "Usage: $0 output-dir"
	exit 1
fi

output=`realpath $1`

cd ${0%/*}

get_control_value()
{
	value=`grep $1 CONTROL | awk '{ print $2 }'`
	echo "${value//[[:space:]]/}"
}

name=$(get_control_value Source)
version=$(get_control_value Version)

#version=`grep Version CONTROL | awk '{ print $2 }'`
#version=${version//[[:space:]]/}

path=$output/$name
filename=$name-$version
zip=$filename.zip

source=$path/$filename

rm -rf $path
mkdir -p $source
cp -r ../*.c $source
cp -r ../*.h $source
cp ../README $source
cp ../CMakeLists.txt $source
cp -r ../Android $source

cp CONTROL $path
sed "s/VERSION/$version/g" portfile.cmake > $path/portfile.cmake

cd $path
zip -9 -r $zip $filename
rm -rf $filename

sha512=`sha512sum $zip | awk '{ print $1 }'`

sed -i "s/CHECKSUM/$sha512/g" portfile.cmake
