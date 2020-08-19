from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired
from wtforms import SubmitField


class FastQForm(FlaskForm):
    forward_file = FileField(validators=[FileRequired()])
    reverse_file = FileField(validators=[FileRequired()])
    submit = SubmitField('Upload')

class TSVForm(FlaskForm):
    tsv_file = FileField(validators=[FileRequired()])
    submit = SubmitField('Upload')
    
class XMLForm(FlaskForm):
    xml_file = FileField(validators=[FileRequired()])
    submit = SubmitField('Upload')
