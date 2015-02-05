#!/bin/bash
#
# Linux utility to switch the default array_tools directory
# to a user-specified one from the default one
#
# Usage: 
#   
#
# Anthony Ho, ahho@stanford.edu, 1/14/2015
# Last update 1/14/2015
#
# WORK IN PROGRESS


# Define usage function
usage() {
    cat <<EOF

$0 

EOF
}

# Define path removal functions
path_remove ()  { 
    PATH2=`echo -n $PATH | awk -v RS=: -v ORS=: '$0 != "'$1'"' | sed 's/:$//'`; 
    echo $PATH2;
}

matlabpath_remove ()  { 
    MATLABPATH2=`echo -n $MATLABPATH | awk -v RS=: -v ORS=: '$0 != "'$1'"' | sed 's/:$//'`; 
    echo $MATLABPATH2
}

# Define function to check and make symbolic links
check_make_sym_links () {
    if [ ! -e $2 ]
    then 
	ln -s $1 $2
    fi
}


# Define paths and directories
central_repo_path="/usr/local/lib/array_tools"
CPlibs="CPlibs"
CPscripts="CPscripts"
sym_bin_path="$HOME/.array_tools_path"

num_arguments=1


# Check number of arguments
if [ $# -gt $num_arguments ]
then

    echo -e "Too many arguments!\n"
    usage
    exit 1

# Switch to new repo if an argument is given 
elif [ $# -eq 1 ]
then

    if [ ! -d "$1" ]
    then
	echo -e "directory given doesn't exist!\n"
	usage
	exit 1
    fi

    new_repo_path=$(readlink -f $1)

    # Creating a hidden directory to store the symbolic links pointing to the new repo
    mkdir -p $sym_bin_path
    
    # Making symbolic links pointing to the new repo
    check_make_sym_links ${new_repo_path}/${CPscripts} ${sym_bin_path}/${CPscripts}
    check_make_sym_links ${new_repo_path}/${CPlibs} ${sym_bin_path}/${CPlibs}
    check_make_sym_links ${new_repo_path}/${CPscripts}/quantifyAllTiles.py ${sym_bin_path}/quantifyAllTiles
    #ln -s 
    
    # Add to path and matlabpath
    export -n PATH="$sym_bin_path:$PATH"
    export -n MATLABPATH="${sym_bin_path}/${CPscripts}:$MATLABPATH"
    export -n MATLABPATH="${sym_bin_path}/${CPlibs}:$MATLABPATH"

# Switch back to default repo if no argument is given
elif [ $# -eq 0 ]
then

    new_repo_path=$(readlink -f $central_repo_path)

    # Remove hidden directory that stores the symbolic links pointing to the new repo
    if [ -d "$sym_bin_path" ]
    then
	rm -r "$sym_bin_path"
    fi

    # Remove to path and matlabpath
    path_remove $sym_bin_path
    matlabpath_remove $sym_bin_path/${CPscripts}
    matlabpath_remove $sym_bin_path/${CPlibs}

fi
