# Developers Guide
While it should work out-of-the box without modification, Arsenal Gear is an
open-source project, and we welcome contributions from the community. This guide
provides information for developers who want to contribute to the project.

## Recommended Development Environment
You are of course free to use whatever set of tools you prefer, but we have
found that the following tools make a fairly low-friction development
environment.  The tools in this section are not completely essential for
developing `arsenal_gear`, but are _strongly_ recommended, and will make your
life and ours easier.  These tools can all be easily installed and set up
automatically on your local machine by running:
```
source setup_dev
```
from the project root directory, and then using the virtual environment that
is built in the project root's `.venv` directory for testing, linting, etc.

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

## Linting and Code Style
We use [pylint](https://pylint.pycqa.org/en/stable/index.html) as our linter.
You can manually run it with:
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
We also use the [black](https://black.readthedocs.io/en/stable/) code formatter
which will make the linter happy at all times provided you run it prior to
running the linter.  You can run black manually with:
```
black arsenal_gear
```

## Automated Tests
We use [pytest](https://pytest.org/) as our testing framework. The testing suite
is located in `tests/`.  You can run them with `pytest` from the project root
directory.

It is good practice to try and get tests that cover as large a fraction of your
code as possible.  One easy way to check your test coverage is with the
`pytest-cov` plugin.  You can install it with `pip install pytest-cov`, and then
run your tests with coverage reporting with:
```
pytest --cov=arsenal_gear tests/
```
This will generate a report showing you the test coverage of each file in the
`arsenal_gear` package using the tests in the `tests/` directory.

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

## Citations
Arsenal Gear wraps up a lot of other people's work: IMFs, yield tables,
isochrones, and stellar evolution models all come from published papers whose
authors deserve to be cited.  Rather than leaving users to work out which
papers a given calculation actually depended on, the code carries that
information itself, and can emit a bibtex bibliography for any entry point.

### Citing a paper in your code
Attach citations with the `cite` decorator, giving it one or more identifiers.
It works on classes and on individual functions or methods:
```python
from arsenal_gear.utils import cite

@cite("10.1086/145971")
class Salpeter(IMF):
    ...
```
Two kinds of identifier are accepted, and the code works out which is which
from the shape of the string, so you never have to say:

- a **DOI**, e.g. `"10.1086/145971"`;
- a **NASA ADS bibcode**, e.g. `"1996A&A...315..105R"`, for papers old enough
  to predate DOIs.  Write the bibcode literally, with a plain `&` and no URL
  escaping.

Pass several identifiers if a class rests on more than one paper, e.g.
`@cite("10.1093/mnras/stz2158", "10.3390/universe7020025")`.  Stacking the
decorator works too, and accumulates rather than replaces.

If a paper is in neither form, it can't be cited here.  Anything NASA ADS
doesn't index isn't really part of the published literature, so there is
deliberately no way to paste in bibtex by hand.

### Regenerating the bibliography
The bibtex itself lives in `arsenal_gear/utils/citations.bib`, which is
generated, not hand-written.  **After adding or changing a `@cite`, regenerate
it and commit the result alongside your code:**
```
python -m arsenal_gear.utils.refresh_citations
```
This imports every `arsenal_gear` module (which is what registers the
citations), looks up whatever the bibliography is missing, drops entries for
citations that no longer exist, and rewrites the file.  Pass `--force` to
re-fetch everything rather than only the missing entries.

The lookup goes to NASA ADS, which has no anonymous API, so you will need a
free API token: sign in at [ADS](https://ui.adsabs.harvard.edu/), generate a
token under Account → Settings → API Token, and put it in your environment as
`ADS_DEV_KEY` (a `.envrc` is a good home for it if you use direnv).  This is a
maintainer-only requirement — because the generated bibliography ships with the
package, users never need a token or a network connection.

The refresher is also where identifiers are checked.  Nothing validates them at
import time, so a typo'd DOI will sit quietly in the source until you run the
refresh, which reports it and exits non-zero.  That's still before anything can
reach a user, since a citation can't resolve until it's in the bibliography.

### Getting a bibliography out
`gather_bibtex` takes a function, method, or class and returns bibtex for
everything it cites:
```python
from arsenal_gear.utils import gather_bibtex
from arsenal_gear.formation.dist_funcs.imf import Salpeter

print(gather_bibtex(Salpeter))
```
```
@ARTICLE{1955ApJ...121..161S,
       author = {{Salpeter}, Edwin E.},
        title = "{The Luminosity Function and Stellar Evolution.}",
      journal = {\apj},
         year = 1955,
...
```
It finds far more than the citations sitting directly on the target.  Starting
from the target it walks the call graph — the names referenced in a function's
body, the methods a class defines, and everything in its MRO — so citations
attached to a base class, to an inherited method, or to a class only reached
indirectly all get collected.  A citation the bibliography is missing becomes a
`% Could not resolve ...` comment rather than an exception, so one stale
reference doesn't cost you the rest of the bibliography.

### How it fits together
The system is split across two modules, and the split is worth preserving:

- `arsenal_gear/utils/citations.py` is what the package uses at runtime.  It
  holds the `cite` decorator, the `REGISTRY` of every identifier seen, the
  call-graph walk in `find_citations`, and the bibliography reader.  It performs
  **no network access** and never writes `citations.bib` — it only reads it.
- `arsenal_gear/utils/refresh_citations.py` is the maintainer tool.  It is the
  only part of the package that talks to NASA ADS, checks identifier shapes, and
  writes the bibliography.

`citations.bib` is an ordinary `.bib` file — you can point LaTeX straight at it
— with a `%%ID <identifier>` comment before each entry recording which citation
it answers.  Bibtex ignores text outside an entry, so those markers cost
nothing.  Note that ADS exports use journal macros (`\apj`, `\mnras`), which
need `aastex` or a similar package to compile.

## GitHub Actions
We use GitHub Actions to run our automated tests on every pull request and push
to the main branch. Currently, we have three workflows:

- `documentation.yml`: This workflow builds the sphinx documentation and deploys it to GitHub Pages.  See [here](docs.md) for more information.
- `pylint.yml`: This workflow runs `pylint` to check for any errors that can be caught by static analysis.
- `pytest.yml`: This workflow runs our unit tests using `pytest`.

## Extra Handy Tools
In addition to the key tools that will be setup when you `source setup_dev`,
there are some useful tools that many of the `arsenal_gear` developers use.

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
