# lpce/tests/conftest.py
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--config-name",
        action="store",
        default="config",
        help="Name of the configuration file",
    )


@pytest.fixture
def config_name(request):
    return request.config.getoption("--config-name")
