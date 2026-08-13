from pcatrfqtl.logging.logger import get_logger


def test_logger_creation():
    logger = get_logger("test")

    assert logger.name == "test"
