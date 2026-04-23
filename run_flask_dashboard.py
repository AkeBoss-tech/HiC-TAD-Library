#!/usr/bin/env python3
from __future__ import annotations

from dotenv import load_dotenv

from flask_dashboard import create_app


load_dotenv()
app = create_app()


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5001, debug=True)
