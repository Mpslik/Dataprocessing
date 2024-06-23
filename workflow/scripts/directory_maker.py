# setup_directories.py
import os
import sys
import yaml


if __name__ == "__main__":
    # Expect the config file path as the first argument
    config_path = sys.argv[1]
    create_directories(config_path)
