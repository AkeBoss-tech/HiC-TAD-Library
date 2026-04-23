from __future__ import annotations

import os
from pathlib import Path

from flask import Flask


def create_app() -> Flask:
    base_dir = Path(__file__).resolve().parent
    app = Flask(
        __name__,
        template_folder=str(base_dir / "templates"),
        static_folder=str(base_dir / "static"),
    )
    app.config["SECRET_KEY"] = os.environ.get("FLASK_SECRET_KEY", "hic-tad-library-dev-secret")
    app.config["PROJECT_ROOT"] = str(base_dir.parent)

    from flask_dashboard.routes import bp

    app.register_blueprint(bp)
    return app
