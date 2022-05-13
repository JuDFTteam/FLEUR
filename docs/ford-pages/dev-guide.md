title: FLEUR developers guidelines

FLEUR software development guide
===================

While the development effort for FLEUR is mainly hosted at [the Institute Quantum Theory of Materials @Forschungszentrum Jülich Germany](https://www.fz-juelich.de/pgi/pgi-1/EN), we welcome contributions to the development from anyone. Please note, that FLEUR is open-source code using the MIT license, so all contributions should adhere to this development model.


In the development of FLEUR we aim at a collaborative, sustainable and flexible development process. The following guidelines are supposed to help in achieving this goal. They are not meant as absolute rules but sketch the general direction. Deviations of the development from this guide can be acceptable if there are good reasons to do so.

## General guidelines

* We reuse code, if new/improved functionality is implemented, instead of developing something new, we try to extend and improve existing code whenever possible. We constantly refactor code to adjust it to new challenges. Ultimately, this should lead to modular code in which functionality is provided that is understandable, reusable and verified and that should be re-useable in sitations not considered explicitely while developing the original code. Try to avoid creating redundant code.

* All new functionality should be accompanied by
   - a test of the functionality. If needed consider also different usage situations, e.g. different levels of parallelism.
   - a documentation that ensures its usability.

* While these guidelines should be taken into account for new code, we have a significant amount of old code. If you touch this, please try to carefully push it into the direction outlined here.

* We have a FLEUR meeting once every two weeks in which we discuss important changes. Please contact d.wortmann@fz-juelich.de for information if needed.

 
## How to write code 

### Modern portable Fortran
  - Encapsulate code into modules, use explicit interfaces to routines where available.
  - Do not use explicit shape dummy arguments without good reason like constant extent.
  - Use data types to collect variables and functions.
      - The collection of variables in data types is done to avoid the creation of subroutine signatures with many arguments and to reduce the need to modify a lot of code if a certain variable is needed in multiple subroutines.
      - The inclusion of procedures and functions into the types can help with modularization and data encapsulation. Functionality that is naturally tied to a specific data type should be provided as a type-bound procedure.
   - Use optional arguments when appropriate to extend the functionality of existing routines.
   
### Code maintainability
  - Have an eye on keeping the required know-how to understand a piece. Try to have simple and local code.
      - Try to avoid writing very large subroutines(probably >1000 lines is too long).
      - Avoid global variables (and also module variables). 
      - Try to write code such that it is as self-explanatory as possible. This includes using meaningful variable names.
  - Aim at a single, clearly defined purpose with each subroutine. Do not sprinkle functionality, e.g. it is not good to just add a few lines of code to a high level routine not at all related to the functionality you are dealing with just because it is convenient.
  - Complex subroutines benefit from comments at the top explaining what the subroutine does and how this is done. 
  - Use the existing functionality, neither re-invent the wheel, nor forget to take advantage of already available features. E.g.:
      - Do not use a simple "STOP" but the judft_error routines as these also cover the case in which multiple MPI processes are used and only a single process triggers an error.
      - Add timers with timerstart/timerstop at appropriate places.
   
### Parallel programming and performance
  - Think about functionality and correctness before focussing on performance.
  - Reuse code already tuned for performance.
     - Use the t_mat type. 
     - Use the MPI wrappers Fleur provides. We do a lot of extra checks there.
  - Nowadays performance is achieved by data-locality, reuse of data and simple data access patterns. 
     - Avoid complex mapping arrays and similar non-local data access.
     - Quite often it is better to improve the performance of the general case instead of having many special code paths which simply try to optimize calculations for special situations. The focus on the general case also often improves maintainabilty and readability of the code.
     - Use optimized math libraries. If possible formulate complex subroutines such that they can be seen as a couple of BLAS calls. Simple subroutines consisting of BLAS calls are often faster than very sophisticated code with many of optimizations. Single High-level BLAS calls are faster than multiple low-level BLAS calls. The same statements also hold for FFT libraries.
  - We currently use three different parallel programming paradigms: MPI, OpenMP, OpenACC.

### Structured Input and Output
  - Do not create 'special' additional files for IO, for switches, for data used in special cases.
  - Use XML IO for stuff that might be processed automatically.
      - Introduce new input parameters and tags as optional in the XML Schema definition of the input file. This simplifies reusage of old input files.
  - Use HDF5 for large binary output.
  - We have have only a single 'out' file as a 'human readable' output and an 'out.xml' file for automatic postprocessing.
  - Try hard to keep IO compatible to previous versions.
 
## Testing
  - Use regression tests, e.g. we have tests to ensure that FLEUR still produces the same results as in older versions. 
  - Take a look at [the docu of the testing system](https://iffgit.fz-juelich.de/fleur/fleur/-/wikis/Testing/Pytest-test-system)
  - Aim for fast-running tests, i.e < 10 sec
      - For this you can depend on the output of any other test already run, or continue from a later point in the process, i.e continue from some density file at the start
      - Use only few k points, reduce cutoff parameters, use small unit cells.

## Documentation
There are several places where documentation can and should be placed.

  - In code documentation: required to ensure that you and other developers can understand your code.
  - [User guide documentation](https://www.flapw.de) (webpage and downloadable files). There is a different [git repository for that](https://iffgit.fz-juelich.de/fleur/www.flapw.de/)
  - [Gitlab wiki](https://iffgit.fz-juelich.de/fleur/fleur/-/wikis/home): Here you can place more in depth information and documentation aimed at developers focussing on specific technical things.

## Usage of the Git repository

   - [Git repository](https://iffgit.fz-juelich.de/fleur/fleur):
      - Small changes should go into the develop branch directly.
      - Alternatively, branches with a short lifetime (until merged into develop again) can be used if additional testing seems appropriate.
      - Large developments that tend to break the develop branch should be done in separate branches. In this case extra care is required that developments do not diverge, e.g. by merging the development branch into your branch regularly.
      - To make tracking of new bugs easy keep single commits small, cover a single purpose with each commit.
      - Use meaningful commit messages. For example see [](https://commit.style/), or [](https://udacity.github.io/git-styleguide/) (we do not use type: in commit messages). As our commits should be small often subjects are sufficient messages.
 

## Code style, specific conventions and hints

   - Use clear function names, that everyone understand (good: `setupHamiltonianOffdiag` bad: `shod`).
   - Use a 3 space indent. No tab indentation.
   - OpenMP
       - Avoid small OpenMP loops but put OpenMP to outer loops when possible. Note that OpenMP can already be beneficial for very basic tasks like copying a large array.
       - If possible use default(none) in OpenMP. This ensures that the code stays maintainable and is not broken by future changes.
   - GPU programming
     -  Use `acc data`-structues rather than `acc enter/exit data` whenever possible.
