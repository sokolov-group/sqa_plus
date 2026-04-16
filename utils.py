import time
import logging
from functools import wraps

def setup_logging(level=logging.INFO):
    logging.basicConfig(
        level=level,
        format="[%(levelname)s] %(name)s: %(message)s",
    )

def log_timing(func):
    func_logger = logging.getLogger(func.__module__)

    @wraps(func)
    def wrapper(*args, **kwargs):
        wall_start = time.perf_counter()
        cpu_start = time.process_time()

        try:
            return func(*args, **kwargs)
        finally:
            cpu_elapsed  = time.process_time() - cpu_start
            wall_elapsed = time.perf_counter() - wall_start

            func_logger.info(
                "%s() wall=%.3fs cpu=%.3fs",
                func.__name__,
                wall_elapsed,
                cpu_elapsed,
            )
    return wrapper

