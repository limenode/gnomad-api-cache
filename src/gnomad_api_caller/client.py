from __future__ import annotations

import requests

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"

def post_gnomad(query: str, variables: dict[str, str]) -> dict:
    """POST a GraphQL query to the gnomAD API and return the JSON response.

    Raises requests.RequestException on transport errors, or KeyError if the
    response has no 'data' field.
    """

    response = requests.post(
        GNOMAD_API_URL,
        json={"query": query, "variables": variables},
        timeout=30,
    )
    response.raise_for_status()
    return response.json()
