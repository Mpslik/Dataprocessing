# setup_directories.py
import os
import sys
import yaml


def create_directories(configuration_path):
    # Load the configuration file
    with open(configuration_path, 'r') as f:
        config = yaml.safe_load(f)

    # Define the list of directories from the config
    directories = [
        config["paths"]["raw_reads"],
        config["paths"]["trimmed_reads"],
        config["paths"]["aligned_reads"],
        config["paths"]["assembled_contigs"],
        config["paths"]["analyzed_metagenome"],
        config["paths"]["visualizations"],
        config["paths"]["coverage"],
        config["paths"]["logs"],
    ]

    # Create directories if they don't exist
    for directory in directories:
        os.makedirs(directory, exist_ok=True)


if __name__ == "__main__":
    # Expect the config file path as the first argument
    config_path = sys.argv[1]
    create_directories(config_path)
