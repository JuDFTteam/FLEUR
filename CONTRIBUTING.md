How to contribute to the development of FLEUR
================

Please consider the following points to contribute to the development of FLEUR:

* Your contributions should agree with our LICENCE. Hence, you have to put your contributions under 
the MIT open-source licence. This essentially means you agree that everyone can use, modify and distribute your code.

* You should use the Git-server at iffgit.fz-juelich.de/fleur/fleur . If you need access, please contact
g.michalicek@fz-juelich.de or d.wortmann@fz-juelich.de.

## Development process

Currently we use the following general policy: if you have small changes or bugfixes these can be committed directly to
the develop branch. Please try to make sure that the CI-pipeline passes.

Larger developments and those breaking the pipelines should create their own branch. But please make sure your code does not diverge. Merge often!


## Coding style
Please read the following comments regarding coding style when you write code to be included in FLEUR.


### General ideas:
In no particular order we try to collect items to consider when writing code for FLEUR


- Instead of 'stop' use calls to judft_error, judft_warn, judft_end
- Do not read and write any files. Files are neither replacements for common-blocks nor a storage for status variables.
Only exceptions:
-- you create nice IO subroutines
-- you write to the typical FLEUR output files

### Code formating
To unify everyones editor settings we agreed on a 3 spaces indent for all of FLEUR.
Hint: (vim)[http://vim.wikia.com/wiki/Converting_tabs_to_spaces] (emacs)[https://www.gnu.org/software/emacs/manual/html_node/efaq/Changing-the-length-of-a-Tab.html]

### Modules
With regard to Fortran modules please follow these ideas:
- modules are named with a prefix 'm_'
- in general each file should contain exactly one module. The names of the files and the modules should correspond.
- the module should start with a 'implicit none' and a 'private' statement. The private statement is particularly important if other modules are used to make sure we do not 'use' from multiple sources. Exceptions are modules collecting use statements with no further code.


### Passing arrays to functions
Fortran offers severals ways to pass arrays to functions. We recommend using either shape-assumed arrays:

```
real, intent(in) :: x(:,:)
```

or allocatable arrays

```
real, intent(in), allocatable :: x(:,:)
```

Due to fortrans legacy you can also often find statements like this(not recommended):

```
real, intent(in) :: x(n_x,n_y)
```

This is not recommended, because it can be (and often is) used to change the rank and size of the array during the call. Consider this small example and try to predict the output:
```
program main 
   implicit none
   real :: x(3,4)
   integer :: i
   
   x = reshape([(i,i=1,12)], [3,4])
   call f(x)
contains
   subroutine f(x)
      implicit none
      real, intent(in) :: x(10)
      write (*,*) x
   end subroutine f
end program main
```

This style is unintuitive and confusing to everyone who tries to understand your code.