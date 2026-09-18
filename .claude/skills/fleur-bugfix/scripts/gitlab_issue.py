#!/usr/bin/env python3
"""Minimal GitLab helper for the FLEUR bugfix workflow.

Talks to the FLEUR project on iffgit.fz-juelich.de (fleur/fleur, project id 19)
using only the Python standard library -- no `requests`, no `glab` needed.

Reading is public.  Creating issues and comments needs a personal access token
with the `api` scope, taken from the first of these that is set:

    FLEUR_GITLAB_TOKEN, GITLAB_TOKEN, GITLAB_API_TOKEN, CI_JOB_TOKEN

Exit codes: 0 ok, 2 usage/API error, 3 no token available (caller should fall
back to letting the user create the issue by hand).

Usage:
    gitlab_issue.py create  --title T (--body B | --body-file F) [--label L]... [--dry-run]
    gitlab_issue.py comment IID       (--body B | --body-file F)                [--dry-run]
    gitlab_issue.py show    IID
    gitlab_issue.py list    [--search S] [--state opened|closed|all]
"""

import argparse
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

HOST = os.environ.get("FLEUR_GITLAB_URL", "https://iffgit.fz-juelich.de").rstrip("/")
PROJECT = os.environ.get("FLEUR_GITLAB_PROJECT", "fleur%2Ffleur")
TOKEN_VARS = ("FLEUR_GITLAB_TOKEN", "GITLAB_TOKEN", "GITLAB_API_TOKEN", "CI_JOB_TOKEN")
API = f"{HOST}/api/v4/projects/{PROJECT}"


def get_token(required):
    for var in TOKEN_VARS:
        tok = os.environ.get(var)
        if tok:
            return tok
    if not required:
        return None
    sys.stderr.write(
        "No GitLab token found (looked at: %s).\n"
        "Create one at %s/-/user_settings/personal_access_tokens with the 'api'\n"
        "scope and export it as FLEUR_GITLAB_TOKEN, or create the issue by hand.\n"
        % (", ".join(TOKEN_VARS), HOST)
    )
    sys.exit(3)


def api(method, path, payload=None, need_token=True):
    token = get_token(required=need_token)
    url = API + path
    data = None
    headers = {"Accept": "application/json"}
    if token:
        headers["PRIVATE-TOKEN"] = token
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", "replace")[:500]
        sys.stderr.write(f"GitLab API {method} {url} failed: {exc.code} {exc.reason}\n{body}\n")
        sys.exit(2)
    except urllib.error.URLError as exc:
        sys.stderr.write(f"Cannot reach {url}: {exc.reason}\n")
        sys.exit(2)


def read_body(args):
    if args.body is not None:
        return args.body
    if args.body_file == "-":
        return sys.stdin.read()
    with open(args.body_file, encoding="utf-8") as fh:
        return fh.read()


def add_body_args(sub):
    grp = sub.add_mutually_exclusive_group(required=True)
    grp.add_argument("--body", help="text (markdown)")
    grp.add_argument("--body-file", help="file holding the text, or - for stdin")
    sub.add_argument("--dry-run", action="store_true", help="print the payload, send nothing")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subs = parser.add_subparsers(dest="cmd", required=True)

    create = subs.add_parser("create", help="open a new issue")
    create.add_argument("--title", required=True)
    create.add_argument("--label", action="append", default=[], help="repeatable, e.g. --label Bug")
    add_body_args(create)

    comment = subs.add_parser("comment", help="add a note to an issue")
    comment.add_argument("iid", type=int)
    add_body_args(comment)

    show = subs.add_parser("show", help="print title/state/url of an issue")
    show.add_argument("iid", type=int)

    lst = subs.add_parser("list", help="list issues")
    lst.add_argument("--search")
    lst.add_argument("--state", default="opened", choices=["opened", "closed", "all"])

    args = parser.parse_args()

    if args.cmd == "create":
        payload = {"title": args.title, "description": read_body(args)}
        if args.label:
            payload["labels"] = ",".join(args.label)
        if args.dry_run:
            print(json.dumps(payload, indent=2))
            return
        issue = api("POST", "/issues", payload)
        print(f"#{issue['iid']} {issue['web_url']}")

    elif args.cmd == "comment":
        payload = {"body": read_body(args)}
        if args.dry_run:
            print(json.dumps(payload, indent=2))
            return
        note = api("POST", f"/issues/{args.iid}/notes", payload)
        print(f"comment {note['id']} added to #{args.iid}: {HOST}/fleur/fleur/-/issues/{args.iid}#note_{note['id']}")

    elif args.cmd == "show":
        issue = api("GET", f"/issues/{args.iid}", need_token=False)
        print(f"#{issue['iid']} [{issue['state']}] {issue['title']}\n{issue['web_url']}\n")
        print(issue.get("description") or "(no description)")

    elif args.cmd == "list":
        query = {"state": args.state, "per_page": "20"}
        if args.search:
            query["search"] = args.search
        issues = api("GET", "/issues?" + urllib.parse.urlencode(query), need_token=False)
        for issue in issues:
            print(f"#{issue['iid']} [{issue['state']}] {issue['title']}")


if __name__ == "__main__":
    main()
