This directory is for unit tests for the Python portions of GeoClaw.  You can
run all of these by running `pytest tests` from `$CLAW/geoclaw` or individually
by 
```
pytest test_[TEST_NAME].py
```
You can also run all of the python tests by using the marker `python` with
```
pytest -m python
```
again in `$CLAW/geoclaw`.

See more details on running and updating tests or debugging issues in the
[Clawpack documentation](https://www.clawpack.org/testing.html).
