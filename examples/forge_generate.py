import os
import sys

from examples.local_generate import main
from esm.sdk import client

if __name__ == "__main__":
    if not os.environ.get("ESM_API_KEY", ""):
        print("Please export your Forge API key as ESM_API_KEY environment variable.")
        sys.exit(1)

    # Run Forge.
    main(client())
