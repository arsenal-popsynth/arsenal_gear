# Developers Guide
While it should work out-of-the box without modification, Arsenal Gear is an
open-source project, and we welcome contributions from the community. This guide
provides information for developers who want to contribute to the project.

## Recommended Development Environment
You are of course free to use whatever set of tools you prefer, but we have
found that the following tools make a fairly low-friction development
environment.

### Virtual Environment
I typically just use the standard python `venv` module to manage my virtual
environments.  I usually store this in `.venv`, as it keeps things out of the
way and tied to the project.  To set it up, you can run:
```
python -m venv .venv
```
And then activate it with:
```
source .venv/bin/activate
```
To install Arsenal Gear and its dependencies, you can then run:
```
pip install -e .[all]
```

The `-e` flag tells pip to install the package in "editable" mode, which just
symlinks the source to the installed package so you don't need to re-run `pip
install` every time you change the source code.

### direnv
Sometimes you need to do more than just set up a python virtual environment, but
also load or change environment variables.  If you don't want to do this
globally, [direnv](https://direnv.net/) is a great tool to manage
directory-specific environment variables.  Once installed, you can create a
`.envrc` that contains anything you would normally put in your `.bashrc` or
`.profile`, and direnv will automatically load it when you `cd` into the
directory.

### VS Code
I'm as opinionated as the next person about my editor of choice, but even an old
vim freak like me has to admit that [VS Code](https://code.visualstudio.com/)
has some handy features.  If you decide to use it, the following extensions can
make your development experience better:
- _direnv_ (`mkhl.direnv`): This lets VS code automatically load your directory-specific environment.
- _GitHub_ Actions (`github.vscode-github-actions`): Monitor the result of GitHub Actions workflows.
- _GitHub_ Copilot (`github.copilot`): LLMs are controversial, but as glorified autocomplete and boilerplate generators they work quite well.
- _Pylance_ (`ms-python.vscode-pylance`): Python language server with type checking and code completion.
- _Python_ (`ms-python.python`): Python support for VS Code.
- _Pylint_ (`ms-python.pylint`): Linting support for Python using pylint.
- _Python Debugger_ (`ms-python.debugpy`): Debugging support for Python.
- _Python Environments_ (`ms-python.vscode-python-envs`): Automatically manage and activate your virtual environments.
- _Vim_ (`vscodevim.vim`): Essential if you have incurable vim brain.

## Linting and Code Style
We use `pylint` as our linter.  You can manually run it with:
```
pylint arsenal_gear/
```
By default, we have it configured to only show `FATAL`, `ERROR`, and `WARNING`
messages (with `exec-used` excluded).  We allow the use of `exec()` because
Arsenal Gear's parameter file is itself a python script.  The security concern
with `exec()` isn't a concern for us because Arsenal Gear is scientific software
that isn't going to be exposed to untrusted input.  If you want to run `pylint`
with all messages enabled (including some of the more annoying convention
messages), you can use:
```
pylint --enable=all arsenal_gear/
```

## Automated Tests
We use `pytest` as our testing framework. The testing suite is located in `tests/`.  You can
run them with `pytest` from the project root directory.

It is good practice to try and get tests that cover as large a fraction of your
code as possible.  One easy way to check your test coverage is with the `pytest-cov` plugin.
You can install it with `pip install pytest-cov`, and then run your tests with coverage reporting with:
```
pytest --cov=arsenal_gear tests/
```
This will generate a report showing you the test coverage of each file in the `arsenal_gear` package using
the tests in the `tests/` directory.

### Pre-Commit Hooks
It's nice to have some tests and linting done automatically before you commit
your changes (this avoids some of the "million little commits to fix a typo"
problem).  We use git's pre-commit hooks to do this.  To set them up, you can run:

```
pip install pre-commit && pre-commit install
```
from within the project directory.  The pre-commit hooks are defined in the
`.pre-commit-config.yaml` file, and will then run every time you run `git
commit`.  Most of the hooks will automatically apply themselves, so if your
commit appears to fail, you can just re-run `git commit` and it will usually
just work.

## GitHub Actions
We use GitHub Actions to run our automated tests on every pull request and push
to the main branch. Currently, we have three workflows:

- `documentation.yml`: This workflow builds the sphinx documentation and deploys it to GitHub Pages.  See [here](docs.md) for more information.
- `pylint.yml`: This workflow runs `pylint` to check for any errors that can be caught by static analysis.
- `pytest.yml`: This workflow runs our unit tests using `pytest`.

### Local Actions
While all our tests can be run locally, you may want to run the tests and/or
linter exactly as it will be run by GitHub.  There is a very cool tool developed
for this called  [act](https://github.com/nektos/act).  `act` uses a Docker
container for a GitHub Actions runner to run your workflows locally.

You will need to set up Docker on your local machine, so follow the instructions
for your operating system [here](https://docs.docker.com/get-docker/).  Docker
is also in many Linux distributions' package managers, so you may be able to
install it with a package manager like `apt` or `dnf` (or load it as a module on
some HPC systems).   Once Docker is set up, you can install `act` following the
instructions [here](https://nektosact.com/installation/index.html).

Running act is straightforward.  To simulate a push, you just run

```
act push
```

and `act` will pull the docker container, spin it up, and run your workflows.
