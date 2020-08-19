from flask import Flask
from flask_bootstrap import Bootstrap
from flask_mail import Mail

from config import Config


mail = Mail()
bootstrap = Bootstrap()


def create_app(config_class=Config):
    app = Flask(__name__)
    app.config.from_object(config_class)

    mail.init_app(app)
    bootstrap.init_app(app)

    from app.api import bp as api_bp
    app.register_blueprint(api_bp, url_prefix='/api')

    from app.web import bp as web_bp
    app.register_blueprint(web_bp)

    return app
