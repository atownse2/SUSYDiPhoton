import os

import git

top_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

repo = git.Repo(top_dir)
git_version = repo.head.object.hexsha[:7]

def print_diff(repo):
    import subprocess
    diff_index = repo.index.diff(None)
    for diff in diff_index.iter_change_type('M'):  # 'M' for modified files
        if diff.a_path.endswith('.ipynb'): continue
        print(subprocess.check_output(['git', 'diff', '--', diff.a_path]).decode())

def check_git():
    repo = git.Repo(top_dir)

    if not repo.is_dirty(): return

    print("Repository has uncommitted changes")
    print_diff(repo)

    _ = input("Would you like to commit these changes now? y/n: ")
    if _ in ['n', 'N']: return
    repo.git.add(update=True)
    message = input("Enter commit message: ")
    repo.index.commit(message)

    print(f"Tags are used to identify versions with distinct outputs, the previous tag is {get_last_tag()}.")
    _ = input(f"Would you like to add a tag to this version? y/n:")
    if _ in ['n', 'N']: return
    tag = input("Enter tag name: ")
    repo.create_tag(tag)

def get_last_tag():
    """Returns the tag associated with the most recent commit"""
    repo = git.Repo(top_dir)
    return repo.tags[-1].name