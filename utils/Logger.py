import logging 

def logger_config(log_path, log_name):

    logger = logging.getLogger(log_name)
    logger.setLevel(level= logging.DEBUG)

    handler = logging.FileHandler(log_path,encoding= 'UTF-8')
    handler.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    handler.setFormatter(formatter)

    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)

    logger.addHandler(handler)
    logger.addHandler(console)

    return logger

if __name__ == "__main__":

    logger = logger_config(log_path = 'log.txt', log_name="啊哈哈哈")

    logger.info("info")
    logger.error("error")
    logger.debug("debug")
    logger.warning("warning")

    logger.info("-"*50)
    logger.info("This is where the analyszer started!")
    logger.info("------------------------------------------")
    logger.info("Making work dir and copy configuration file to it....")
    print("for logging!!")