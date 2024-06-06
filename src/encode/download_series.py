import argparse
import os.path
import sys

import requests


def download_carts(cart_name, save_dir):
    if not os.path.exists(save_dir):
        raise ValueError(f'Directory "{save_dir}" doesn\'t exist. Please create it first')

    cart = requests.get(f"https://www.encodeproject.org/carts/{cart_name}/?format=json", timeout=5)
    if not cart.ok:
        raise RuntimeError(f"Could not retrieve cart {cart_name}: {cart.status_code} {cart.raise_for_status}")

    cart = cart.json()

    for element in cart["elements"]:
        if not element.startswith("/functional-characterization-series"):
            continue

        series_id = element[36:-1]

        # https://www.encodeproject.org/functional-characterization-series/ENCSR336WSZ/?format=json
        series_url = f"https://www.encodeproject.org{element}?format=json"
        response = requests.get(series_url, timeout=3)
        if not response.ok:
            print(f"{series_id}: {response.status_code}")
            sys.exit(1)

        with open(os.path.join(save_dir, f"{series_id}.json"), "wb") as out_file:
            out_file.write(response.content)


def get_args():
    parser = argparse.ArgumentParser(description="Download Characterization Series from an Encode cart")
    parser.add_argument("cart_name", help="The name of the cart to download series from")
    parser.add_argument(
        "-s",
        "--save-dir",
        help="The directory to save the cart files in",
        default="temp_data",
    )
    return parser.parse_args()


def run_cli():
    args = get_args()
    try:
        download_carts(args.cart_name, args.save_dir)
    except ValueError as ve:
        print(ve)
        sys.exit(1)
