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

---

## Additional steps (SSH)

### Generate SSH key

`ssh-keygen -t rsa`

### Add the public key to your Github account

In the upper-right corner of any page, click your profile photo, then click Settings.

![]({{site.baseurl}}/https://help.github.com/assets/images/help/settings/userbar-account-settings.png)

![](https://help.github.com/assets/images/help/settings/userbar-account-settings.png)

In the user settings sidebar, click SSH and GPG keys.

![](https://help.github.com/assets/images/help/settings/settings-sidebar-ssh-keys.png)

Click New SSH key or Add SSH key.

![](https://help.github.com/assets/images/help/settings/ssh-add-ssh-key.png)

Paste your key into the "Key" field.

### Test the connection

`ssh -T git@github.com`

### Clone using SSH

### Direct push without entering the credentials

---
---

## Resources

### Cheat Sheets

[GIT Cheat Sheet](https://rogerdudler.github.io/git-guide/files/git_cheat_sheet.pdf)

[MD Cheat Sheet](https://www.markdownguide.org/cheat-sheet/)

### Markdown Live Editors

[StackEdit](https://stackedit.io)

[HackMD](https://hackmd.io)

### UI Software

[gitKraken](https://www.gitkraken.com/)

