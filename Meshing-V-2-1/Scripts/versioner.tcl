#!/bin/sh
# the next line restarts using tclsh \
    exec tclsh "$0" "$@"

################################################################################
#
# versioner.tcl
#
# adds a lib/Slicer3/Slicer3Version.txt file to the current build
#
# Usage: (no options)
#   versioner 
#
# Initiated - sp - 2009-03-24
#

################################################################################

# for subversion repositories (Sandbox)
if {[info exists ::env(SVN)]} {
    set ::SVN $::env(SVN)
} else {
    set ::SVN svn
}


################################################################################
# build the lib/Slicer3/Slicer3version.txt file
# - determine the location
# - determine the build variables
# - get the svn info
# - write the file
# 

set cwd [pwd]
cd [file dirname [info script]]
cd ..
set ::Slicer3_HOME [pwd]

source $::Slicer3_HOME/slicer_variables.tcl

set ::Slicer3_BUILDDATE [clock format [clock seconds] -format %Y-%m-%d]

set svninfo [split [exec svn info] "\n"]
array set svn ""
foreach line $svninfo {
  foreach {tag value} $line {
    if { $tag == "URL:" } {
      set svn(URL) $value
    }
    if { $tag == "Revision:" } {
      set svn(revision) $value
    }
  }
}
cd $cwd

set fp [open $::Slicer3_BUILD/lib/Slicer3/Slicer3Version.txt "w"]
puts $fp "build $::env(BUILD)"
puts $fp "buildDate $::Slicer3_BUILDDATE"
puts $fp "svnurl $svn(URL)"
puts $fp "svnrevision $svn(revision)"
close $fp
