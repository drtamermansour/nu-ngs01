# GIT version control tutorial

## Creation

1. Register and signin on [Github](https://www.github.com)
2. Create a new repository with the name **temp-repo**

### 3- Clone the repo on your own device

```bash
git clone https://github.com/mr-eyes/temp-repo.git
cd temp-repo/
```

### 4- Create README.md and update its content

`echo "# Welcome to the NGS temp repo" > README.md`

### 5- Sync the changes

```bash
git add README.md
git commit -m "create readme"

# you may need to enter your credentials
git push
```

### 6- Create another file on the web UI

### 7- Pull the changes locally

`git pull`

### 8- Create a new branch and edit the second file

### 9- Compare branches

### 10- Create Pull Request and merge

---

## Quick Team Project

### Group yourselves to 3:5-members groups.

### One of you create a repo and the others fork and clone it.

### Other members modify some file and create a PR to the creator's branch.

### The creator accept one of the previous PRs.

### Try the collaboratos features for direct pushes to the repo.