
##########################
Testing Driven Development
##########################

A typical software development process lifecycle start with the development of unittests, developing code for the unittests, proceeding to integration testing to deal with compatibility of code between two or more different units.  Finally, a functional test is run to ensure that the application produces the type of responses which we expect it to.  

Unit Testing
============
Unit testing is the act of testing a "unit" in an application.  Within this context, a "unit" is often a function or method of a class instance.  This unit is referred to the "unit under test."  The goal of unit testing is to test permutations of the unit-under-test.  

The major benefit of unit testing is that it identifies problems early in the development lifecycle.

Unit testing in native within the Python programming language, and is contained in the module :obj:`unittest`.  A lot has been written on the process of unittesting.

Integration Testing
===================

Integration testing tests two different parts of the application at once,
including the integration of the parts, to determine if they function as 
intended.  This type of testing identifies defects in the interfaces between
disparate parts of the codebase as well as interaction to other parts of the 
program.

Functional Testing
==================

Problems with this Approach
###########################
This software development process is actually too heavy for a one-person development team or for a research group where individuals are mostly working independently.  From a practical standpoint, unittests need to be written when code becomes shared because unittests become a contract between the developer and users because it serves as a reference implmentation of specific methods and clearly defines their behavior.  However, the development of this code

