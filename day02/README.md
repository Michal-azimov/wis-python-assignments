# Day 02 – Issues, Python, and Hello World

## 1. Reviewing classmates' web sites with issues

I looked at the [list of students](https://github.com/Code-Maven/wis-python-course-2026-10/blob/main/STUDENTS.md)
and picked two whose site repositories had no issues yet:

- **Delfina Papinutti**: <https://delfip.github.io>. Issue: ISSUE_LINK_1
- **Ehud Elfer**: <https://ehudelfer.github.io>. Issue: ISSUE_LINK_2

I read each site and opened an issue on its repository with praise and
suggestions for improvement. Then I followed the discussion and responded to
issues opened on my own site,
[Michal-azimov.github.io](https://github.com/Michal-azimov/Michal-azimov.github.io/issues).
When I changed my site because of an issue, I put the issue number in the commit
message (for example `Improve layout, fixes #1`) and closed the issue when I
was done.

## 2. Installing Python

I installed Python 3 on Windows and checked the installation from the terminal:

```
python --version
```

## 3. Hello World in VS Code

I opened this folder in VS Code and created [`hello.py`](hello.py):

```python
print("Hello World")
```

I ran it from the VS Code terminal:

```
python day02/hello.py
```

Output:

```
Hello World
```

## 4. The assignments repository

1. Created a new public repository on GitHub for my assignments.
2. Cloned it to my computer with `git clone`.
3. Added a top-level `README.md` describing the repository.
4. Created the `day01` folder with a `README.md` describing the first assignment.
5. Created the `day02` folder with `hello.py` and this `README.md`.
6. Committed the files and pushed them to GitHub (`git add`, `git commit`, `git push`).
7. Opened the issue "Day 02 by Michal Azimov" on the course repository, with
   links to this repository and to the two issues I opened.
