# Developers Guide
While it should work out-of-the box without modification, Arsenal Gear is an
open-source project, and we welcome contributions from the community. This guide
provides information for developers who want to contribute to the project.

## Recommended Development Environment
You are of course free to use whatever set of tools you prefer, but we have found that the following
tools make a fairly low-friction development environment.

### Virtual Environment
### direnv
Sometimes you need to do more than just set up a python virtual environment, but also load or change
environment variables.  If you don't want to do this globally, [direnv](https://direnv.net/) is a great tool
to manage directory-specific environment variables.  Once installed, you can create a `.envrc` that contains
anything you would normally put in your `.bashrc` or `.profile`, and direnv will automatically load it when you
`cd` into the directory.
### Pre-Commit Hooks
### VS Code
I'm as opinionated as the next person about my editor of choice, but even an old vim freak like me has to admit that
[VS Code](https://code.visualstudio.com/) has some handy features.  If you decide to use it, the following extensions
can make your development experience better:
- direnv (mkhl.direnv): This lets VS code automatically load your directory-specific environment.
- Github Actions (github.vscode-github-actions): Monitor the result of Github Actions workflows.
- Github Copilot (github.copilot): LLMs are controversial, but as glorified autocomplete and boilerplate generators they work quite well.
- Pylance (ms-python.vscode-pylance): Python language server with type checking and code completion.
- Python (ms-python.python): Python support for VS Code.
- Pylint (ms-python.pylint): Linting support for Python using pylint.
- Python Debugger (ms-python.debugpy): Debugging support for Python.
- Python Environments (ms-python.vscode-python-envs): Automatically manage and activate your virtual environments.
- Vim (vscodevim.vim): Essential if you have incurable vim brain.

## Linting and Code Style

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

## Github Actions
We use Github Actions to run our automated tests on every pull request and push to the main branch. Currently, we have three workflows:

- `documentation.yml`: This workflow builds the sphinx documentation and deploys it to GitHub Pages.  See [here](docs.md) for more information.
- `pylint.yml`: This workflow runs `pylint` to check for any errors that can be caught by static analysis.
- `pytest.yml`: This workflow runs our unit tests using `pytest`.

### Local Actions
