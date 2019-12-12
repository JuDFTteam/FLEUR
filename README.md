<div align="center">
<img src="https://www.flapw.de/site/img/fleur.gif"  width="220">


Welcome to the source code of FLEUR
=====================

[Report bug](https://iffgit.fz-juelich.de/fleur/fleur/issues/new?template=Bug.md)
.
[Request feature](https://iffgit.fz-juelich.de/fleur/fleur/issues/new?template=FeatureRequest.md&labels=feature)

[Homepage and Documentation](https://www.flapw.de)
</div>

## Table of contents

- [Using the FLEUR git repository](#fleur-git-repository)
- [Dealing with Bugs and problems](#bugs-and-problems)
- [Installation of FLEUR](#installation-of-FLEUR)
- [Contributing](#contributing)

## FLEUR git repository

The primary git-repository of FLEUR can be found on the [iffgit-Server at FZ-Jülich](https://iffgit.fz-juelich.de/fleur/fleur/).

You can clone the repository by using
```
git clone https://iffgit.fz-juelich.de/fleur/fleur.git
```

If you are a FLEUR developer you should use
```
git clone gitlab@iffgit.fz-juelich.de:fleur/fleur.git
```
to be able to push changes back to the server. If you are not a developer yet but want to contribute, please contact [Gregor](g.michalicek@fz-juelich.de) or [Daniel](d.wortmann@fz-juelich.de).

Please note, that the default branch you will see after cloning the repository is the 'develop' branch. In general you might find 
the following branches on the server.

* develop: this is the default branch with the most up-to-date version of FLEUR. Small changes and developments should be committed 
directly into this branch. When doing so you should try to keep the code operational. It should still compile and the test should run. 
* release: this branch collects the official releases. You cannot commit to this branch and bugfixes should be handled as [described below](#bugs-and-problems).
* stable: this branch contains snapshots of the development branch considered "stable".

In addition several other branches can/will be present corresponding to features currently under development. If you start your own larger development
it can be advisable to create your own branch. In this case you should try to follow changes in 'develop' by frequently merging 'develop' into your branch
and you should create a merge request with 'develop' as soon as you are finished or reached some usefull state in your development.

## Bugs and Problems

You might experience bugs in FLEUR :-).

If you find a bug you should:

A)  [Report this bug by generating an Issue](https://iffgit.fz-juelich.de/fleur/fleur/issues/new?template=Bug.md). Please describe in 
detail the relevant input and what happens. You should consider using 
the bug-template for your issue as this will help you providing us with 
the relevant information.

or/and

B) Provide a bugfix. If the bug is only present in the development branch/ is due
to a new feature under development simply commit your fix to the development branch.
If you are fixing a bug in a release-version, please:
* check out the git release branch: ```git checkout --track origin/release```
* create a bugfix branch: ```git checkout -b bugfix_SOME_NAME_HERE```
* apply your changes, test them and commit them
* push your bugfix branch to the server: ``` git push -u origin bugfix_SOME_NAME_HERE```
* create a merge request on the gitlab to have you bugfix merged with the release branch
* check out the develop branch: ```git checkout develop```
* merge your fix into the develop branch: ```git merge bugfix_SOME_NAME_HERE```


## Installation of FLEUR

To install and use FLEUR, please check the [Documentation](https://www.flapw.de).

## Contributing

FLEUR is an open source code under the [MIT-license](https://iffgit.fz-juelich.de/fleur/fleur/blob/develop/LICENSE).

Your are very welcome to contribute to its development. If you need help or access to the git repository, 
please contact [Gregor](g.michalicek@fz-juelich.de) or [Daniel](d.wortmann@fz-juelich.de).

Please also use the [Wiki](https://iffgit.fz-juelich.de/fleur/fleur/wikis/home) for sharing information relevant for developers. 