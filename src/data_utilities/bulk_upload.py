#!/usr/bin/env python

import argparse
import csv
import re
import time
from dataclasses import dataclass
from typing import Optional, TypeAlias

import requests

JSON_MIME = "application/json"

BLANK_ACCESSION = ["", None]

URLString: TypeAlias = str
CSRFToken: TypeAlias = str


def print_response(response, title=""):
    print("************************")
    print(title)
    print(response.status_code, response.reason)
    print(response.text[:600])
    print("************************")


def is_url(path: str):
    return path.startswith("https://") or path.startswith("http://")


@dataclass(init=False)
class Metadata:
    accession_id: Optional[str]
    compressed_data_url: Optional[str] = None
    compressed_data_file: Optional[str] = None
    experiment_metadata_url: Optional[str] = None
    experiment_metadata_file: Optional[str] = None
    analysis_metadata_url: Optional[str] = None
    analysis_metadata_file: Optional[str] = None

    def __init__(
        self,
        accession_id,
        experiment_path: Optional[str] = None,
        analysis_path: Optional[str] = None,
        compressed_path: Optional[str] = None,
    ):

        self.accession_id = accession_id

        if experiment_path is not None:
            if is_url(experiment_path):
                self.experiment_metadata_url = experiment_path
            else:
                self.experiment_metadata_file = experiment_path

        if analysis_path is not None:
            if is_url(analysis_path):
                self.analysis_metadata_url = analysis_path
            else:
                self.analysis_metadata_file = analysis_path

        if compressed_path is not None:
            if is_url(compressed_path):
                self.compressed_data_url = compressed_path
            else:
                self.compressed_data_file = compressed_path

    def __str__(self):
        experiment = (
            self.experiment_metadata_file if self.experiment_metadata_file is not None else self.experiment_metadata_url
        )
        experiment_string = "" if experiment is None else f" {experiment}"

        analysis = (
            self.analysis_metadata_file if self.analysis_metadata_file is not None else self.analysis_metadata_url
        )
        analysis_string = "" if analysis is None else f" {analysis}"

        compressed = self.compressed_data_file if self.compressed_data_file is not None else self.compressed_data_url
        compressed_string = "" if compressed is None else f" {compressed}"

        return f"{self.accession_id}: {experiment_string} {analysis_string} {compressed_string}"


class UploadSession:
    hostname: URLString
    login_url: URLString
    logout_url: URLString
    upload_url: URLString
    status_url: URLString

    def __init__(self, hostname, username, password):
        self.hostname = hostname
        self.login_url = f"{hostname}/accounts/login/"
        self.logout_url = f"{hostname}/accounts/logout/"
        self.upload_url = f"{hostname}/uploads/"
        self.status_url = f"{hostname}/task_status/task/"
        self.session = None
        self.username = username
        self.password = password
        self.login_csrf = None
        self.upload_csrf = None
        self.logout_csrf = None

    def __enter__(self):
        self.session = requests.Session()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.session is None:
            raise RuntimeError("Please use UploadSession as a Context Manager")

        self.session.close()

    def _retry_request(self, request, check_json=False):
        response = None
        for _ in range(3):
            response = request()

            if response.ok:
                if check_json:
                    try:
                        status_result = response.json()
                        return status_result
                    except:
                        # The log token doesn't last forever!
                        if "Log In" in response.text:
                            self.login()
                else:
                    return response

            if response.status_code == 403:
                self.login()

            time.sleep(3)

        print_response(response, "Request Error")
        raise RuntimeError("Connection Error")

    def _check_response_status(self, response: requests.Response):
        if not response.ok:
            raise RuntimeError(f"Bad Response Code: {response.status_code} ({response.request.method} {response.url})")

    def _get_login_csrf(self):
        if self.session is None:
            raise RuntimeError("Please use UploadSession as a Context Manager")

        if self.login_csrf is not None:
            # There seems to be an issue where a user's session has expired, but they aren't
            # actually logged out. This causes problems when trying to get a new login csrf token
            # So we'll just make sure they are logged out.
            logout_response = self._retry_request(lambda: self.session.get(self.logout_url))
            search_result = re.search(r'name="csrfmiddlewaretoken" value="(.*)"', logout_response.text)
            if search_result is not None:
                self.logout_csrf = search_result.group(1)
            else:
                print_response(logout_response, "Logout CSRF Request Error")

            self.session.post(
                self.logout_url,
                data={
                    "csrfmiddlewaretoken": self.logout_csrf,
                },
                headers={"Referer": self.logout_url},
            )

        response = self._retry_request(lambda: self.session.get(self.login_url))
        search_result = re.search(r'name="csrfmiddlewaretoken" value="(.*)"', response.text)
        if search_result is not None:
            self.login_csrf = search_result.group(1)
        else:
            print_response(response, "Login CSRF Request Error")

    def _get_upload_csrf(self):
        if self.session is None:
            raise RuntimeError("Please use UploadSession as a Context Manager")

        response = self._retry_request(lambda: self.session.get(self.upload_url))
        search_result = re.search(r'name="csrfmiddlewaretoken" value="(.*)"', response.text)
        if search_result is not None:
            self.upload_csrf = search_result.group(1)
        else:
            print_response(response, "Upload CSRF Request Error")

    def login(self):
        if self.session is None:
            raise RuntimeError("Please use UploadSession as a Context Manager")

        self._get_login_csrf()
        login = self.session.post(
            self.login_url,
            data={
                "csrfmiddlewaretoken": self.login_csrf,
                "login": self.username,
                "password": self.password,
            },
            headers={"Referer": self.login_url},
        )
        self._check_response_status(login)

    def upload(self, metadatum: Metadata):
        if self.session is None:
            raise RuntimeError("Please use UploadSession as a Context Manager")

        self._get_upload_csrf()
        print(f"Uploading {metadatum}")
        data = {
            "csrfmiddlewaretoken": self.upload_csrf,
            "accept": JSON_MIME,
            "experiment_accession": metadatum.accession_id,
        }
        files = {}

        #
        # Add the experiment and analysis metadata to the request as appropriate
        #
        if metadatum.experiment_metadata_file is not None:
            if isinstance(metadatum.experiment_metadata_file, str):
                files["experiment_file"] = open(metadatum.experiment_metadata_file, encoding="utf-8")

        if metadatum.experiment_metadata_url is not None:
            if isinstance(metadatum.experiment_metadata_url, str):
                data["experiment_url"] = metadatum.experiment_metadata_url

        if metadatum.analysis_metadata_file is not None:
            if isinstance(metadatum.analysis_metadata_file, str):
                files["analysis_file"] = open(metadatum.analysis_metadata_file, encoding="utf-8")

        if metadatum.analysis_metadata_url is not None:
            if isinstance(metadatum.analysis_metadata_url, str):
                data["analysis_url"] = metadatum.analysis_metadata_url

        if metadatum.compressed_data_file is not None:
            if isinstance(metadatum.compressed_data_file, str):
                files["full_file"] = open(metadatum.compressed_data_file, "rb")

        if metadatum.compressed_data_url is not None:
            if isinstance(metadatum.compressed_data_url, str):
                data["full_url"] = metadatum.compressed_data_url

        #
        # Upload the metadata
        #
        upload_response = self._retry_request(
            lambda: self.session.post(self.upload_url, data=data, files=files, headers={"Referer": self.login_url}),
            check_json=True,
        )

        task_status_id = upload_response["task_status_id"]

        if (file := files.get("experiment_file")) is not None:
            file.close()
        if (file := files.get("analysis_file")) is not None:
            file.close()
        if (file := files.get("full_file")) is not None:
            file.close()

        #
        # Check the status until the upload is finished or errors out
        #
        status_url = self.status_url + task_status_id
        print(f"Checking task status at {status_url}")
        statuses = set()
        while True:
            status_result = self._retry_request(
                lambda: self.session.get(status_url, params={"accept": JSON_MIME}),
                check_json=True,
            )

            status = status_result["status"]

            match status:
                case "W":
                    message = "Waiting to Start"
                case "S":
                    message = "Started"
                case "F":
                    message = "Finished Successfully"
                case "E":
                    message = f"Server Error: {status_result['error_message']}"
                case _:
                    raise ValueError(f'Unknown status "{status}"')

            if status not in statuses:
                print(message)
                statuses.add(status)

            match status:
                case "F" | "E":
                    break
                case _:
                    time.sleep(3)


def full_hostname(hostname: str) -> URLString:
    if hostname[-1] == "/":
        hostname = hostname[:-1]

    if is_url(hostname):
        return hostname

    return f"http://{hostname}"


def next_accession(accession: str) -> str:
    match = re.fullmatch(r"DCPEXPR([\da-fA-F]{10})", accession)

    if match is None:
        raise ValueError(f"Could not match accession {accession}")

    number = int(match.group(1), base=16) + 1

    return f"DCPEXPR{number:010X}"


# TSV of metadata:
#
# Accession ID\texperiment metadata\tanalysis metadata
#
# If an accession id is not set, use the previous accession id + 1
# Skip blank lines
# Skip lines that start with '#'
# experiment and analysis files are optional
#
def read_metadata_list(metadata_file: str) -> list[Metadata]:
    metadata = []
    prev_accession = None
    with open(metadata_file, encoding="utf-8") as metadata_tsv:
        reader = csv.reader(metadata_tsv, delimiter="\t")
        for row in reader:
            # Skip if the row doesn't have the right number of columns
            try:
                accession, experiment_path, analysis_path = row
            except ValueError:
                continue

            # skip if the row is blank, minus the tabs
            if accession == experiment_path == analysis_path == "":
                continue

            # skip if the row is commented out
            if accession.startswith("#"):
                continue

            # If the row doesn't define an accession ID figure out what
            # it should be
            if accession in ["", None] and prev_accession not in BLANK_ACCESSION:
                accession = next_accession(prev_accession)

            prev_accession = accession

            if len(experiment_path) == 0:
                experiment_path = None

            if len(analysis_path) == 0:
                analysis_path = None

            metadata.append(
                Metadata(
                    accession_id=accession,
                    experiment_path=experiment_path,
                    analysis_path=analysis_path,
                )
            )

    return metadata


# TSV of compressed full data sets:
#
# Accession ID\tcompressed file
#
# If an accession id is not set, use the previous accession id + 1
# Skip blank lines
# Skip lines that start with '#'
#
def read_compressed_list(metadata_file: str) -> list[Metadata]:
    metadata = []
    prev_accession = None
    with open(metadata_file, encoding="utf-8") as metadata_tsv:
        reader = csv.reader(metadata_tsv, delimiter="\t")
        for row in reader:
            # Skip if the row doesn't have the right number of columns
            try:
                accession, compressed_path = row
            except ValueError:
                continue

            # skip if the row is blank, minus the tabs
            if accession == compressed_path == "":
                continue

            # skip if the row is commented out
            if accession.startswith("#"):
                continue

            # If the row doesn't define an accession ID figure out what
            # it should be
            if accession in ["", None] and prev_accession not in BLANK_ACCESSION:
                accession = next_accession(prev_accession)

            prev_accession = accession

            print(f"accession: {accession}")

            metadata.append(Metadata(accession_id=accession, compressed_path=compressed_path))

    return metadata


def get_metadata(args: argparse.Namespace) -> list[Metadata]:
    if args.mfile is not None:
        metadata = read_metadata_list(args.mfile)
    elif args.cfile is not None:
        metadata = read_compressed_list(args.cfile)
    else:
        metadata = [
            Metadata(
                accession_id=args.accession,
                experiment_path=args.ef,
                analysis_path=args.af,
            )
        ]
    return metadata


def get_args():
    parser = argparse.ArgumentParser(description="A script for bulk-uploading experiments/screens to the CCGR Portal")
    parser.add_argument("host", help="Portal hostname")
    parser.add_argument("-u", "--username", help="Username", required=True)
    parser.add_argument("-p", "--password", help="Password", required=True)

    bulk = parser.add_argument_group(title="Bulk Experiment/Analysis Metadata Upload")
    bulk.add_argument(
        "-m",
        "--mfile",
        help="A TSV containing locations of the experiment/analysis metadata needed to load the data",
    )

    bulk = parser.add_argument_group(title="Bulk Full Experiment/Analysis Upload")
    bulk.add_argument(
        "-c",
        "--cfile",
        help="A TSV containing locations of full (compressed) sets of experiment/analysis metadata and data",
    )

    single = parser.add_argument_group(title="Single Experiment/Analysis Upload")
    single.add_argument("-a", "--accession", help="The experiment accession ID")
    single.add_argument("--ef", "--experiment-file", help="The experiment metadata")
    single.add_argument("--af", "--analysis-file", help="The analysis metadata")

    ns = parser.parse_args()

    if (
        int(getattr(ns, "mfile", None) is not None)
        + int(getattr(ns, "cfile", None) is not None)
        + int(
            getattr(ns, "accession", None) is not None
            or getattr(ns, "ef", None) is not None
            or getattr(ns, "af", None) is not None
        )
    ) != 1:
        parser.error("Bulk (-m), Compressed (-c) and Single (-a, --ef, --af) Uploads are mututally exclusive")

    return ns


def main(args: argparse.Namespace):
    host = full_hostname(args.host)
    metadata = get_metadata(args)

    for metadatum in metadata:
        with UploadSession(host, args.username, args.password) as session:
            session.login()
            session.upload(metadatum)


def run_cli():
    parsed_args = get_args()
    main(parsed_args)
