#Adding an existing project to GitHub using the command line

Simple steps to add existing project to Github.

## 1. Create a new repository on GitHub.
In Terminal, change the current working directory to your local project.

## 2. Initialize the local directory as a Git repository.

	git init

Create .gitignore file and add rule to ignore data from commit

    touch .gitignore

Add the files in your new local repository. This stages them for the first commit.

	git add .

Commit the files that you've staged in your local repository.

	git commit -m 'First commit'
    git commit -m 'add md files'
    git commit -m 'v1.0 preserve dir structure'
	git commit -m 'v2.0 sunday work'
	git commit -m 'v4.0 update output/plot-table'



Copy remote repository URL field from your GitHub repository, in the right sidebar, copy the remote repository URL.

In Terminal, add the URL for the remote repository where your local repostory will be pushed.

	git remote add origin https://github.com/thehung92/asd_wes.git
	
Sets the new remote:
	
	git remote -v

Push the changes in your local repository to GitHub.
	`git push origin master
	git push`

Pushes the changes in your local repository up to the remote repository you specified as the origin

## remove remote
	`git remote remove origin`
## git revert change

	git reset HEAD~1

## remove files from tracking in terminal if .gitingore don't work

	find . -name .DS_Store -print0 | xargs -0 git rm --ignore-unmatch
	find . -name *.Rproj -print0 | xargs -0 git rm --ignore-unmatch

