#!/usr/bin/env python

import argparse
import os.path
import sys

import requests


def run(cart_name):
    cart = requests.get(f"https://www.encodeproject.org/carts/{cart_name}/?format=json")
    if not cart.ok:
        raise RuntimeError(
            f"Could not retrieve cart {cart_name}: {cart.status_code} {cart.raise_for_status}"
        )

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

        with open(os.path.join("temp_data", f"{series_id}.json"), "wb") as out_file:
            out_file.write(response.content)


def get_args():
    parser = argparse.ArgumentParser(
        description="Download Characterization Series from an Encode cart"
    )
    parser.add_argument(
        "cart-name", help="The name of the cart to download series from"
    )
    return parser.parse_args()


# to format all these files i ran
#  for f in series_data/*.json; do cat $f | jq '.' > $f.fmt; done
#  for f in series_data/*; do mv $f ${f%.*}; done
# Or close to it.

if __name__ == "__main__":
    args = get_args()
    run(getattr(args, "cart-name"))
