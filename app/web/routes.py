import datetime
import os
import re
import uuid
from shutil import rmtree

from flask import render_template, redirect, url_for, session, flash, request
from werkzeug.utils import secure_filename

from app.core.preprocessing.parser import preprocess_fastq_files, initialize_empty_dataframe, preprocess_tsv_file, splitTSV
from app.core.utils.preprocess_utils import allowed_file
from app.core.objects.FastQDataframe import FastQDataframe
from app.web import bp
from app.web.forms import FastQForm, TSVForm, XMLForm

#last_purge = None

#@bp.before_request
#def before_request():
    #global last_purge
    #now = datetime.datetime.now()

    #try:
        #if last_purge:
            #delta = now - last_purge
            #if delta.seconds > 3600:
                #for dir in os.listdir(os.path.abspath("data/")):
                    #if re.search("([a-zA-Z0-9_]+-)+([a-zA-Z0-9_]+)", dir):
                        #for file in os.listdir(os.path.abspath("data/"+dir)):
                            #print(os.path.join(os.path.abspath("data/"+dir), file), "automatically removed")
                            #os.remove(os.path.join(os.path.abspath("data/"+dir), file))

                #last_purge = now
        #else:
            #last_purge = now
    #except:
        #pass

@bp.route('/')
@bp.route('/index')
def index():
    return render_template('index.html')

#Function that will load the preprocess page and collect stats for on the page.
@bp.route('/preprocessing', methods = ['GET','POST'])
def preprocessing():
    session_id = session['id']
    try:
        fastq_df = FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl")
        gotBlast = "fw_full_tax" in fastq_df.get_dataframe().columns
    except FileNotFoundError:
        return render_template('expired.html')
    dictStats = {}
    fwTotal, rvTotal = fastq_df.getTotalReads()
    minLength = fastq_df.getMinLength()
    maxLength = fastq_df.getMaxLength()
    maxForwardLength = fastq_df.getMaxLength("fw")
    maxReverseLength = fastq_df.getMaxLength("rv")
    fwMaxQuality, rvMaxQuality = fastq_df.getMaxQuality()

    #Fill a python dictionary with stats about the input data.
    dictStats.update(fastq_df.get_preprocessOptions())
    dictStats.update({"fw_total": fwTotal, "rv_total": rvTotal, "min_length": minLength, "max_length": maxLength, "max_fw_length": maxForwardLength, "max_rv_length": maxReverseLength, "fw_max_quality": fwMaxQuality,
                      "rv_max_quality": rvMaxQuality
})    
    jsonOptions = fastq_df.get_originalPreprocessOptionsJSON()

    return render_template('preprocessing.html', stats = jsonOptions, stats2 = dictStats, gotBlast = gotBlast)

# Function that provide the buttons to import files and calls functions to process the input
@bp.route('/indexNewSession', methods=['GET', 'POST'])
def indexNewSession():
    form = FastQForm()            # This variable is used to create the import buttons and their functionality

    # Code that is used when the user press the submit button:
    if form.validate_on_submit():

        # Functions to check if the file type is correct:
        fw_file = secure_filename(form.forward_file.data.filename)
        rv_file = secure_filename(form.reverse_file.data.filename)

        # When the files are of the correct filetype use the following code:
        if allowed_file(fw_file) and allowed_file(rv_file):

            try:
                session_id = get_sessionID()
                if os.path.exists('data/'+session_id):
                    if session.get("id"):
                        session.pop("id")
                    session_id = get_sessionID()
            finally:
                # Rename the files:
                renamed_fw_file = 'fw_file.'+fw_file.rsplit('.', 1)[1].lower()
                renamed_rc_file = 'rv_file.'+rv_file.rsplit('.', 1)[1].lower()

                # Create directory and save the files in it:
                if not os.path.exists('data/'+session_id):
                    try:
                        os.makedirs('data/'+session_id)
                        form.forward_file.data.save('data/' + session_id + '/' + renamed_fw_file)
                        form.reverse_file.data.save('data/' + session_id + '/' + renamed_rc_file)

                        # Create the FastQDataFrame object with the files and export it as a pickle:
                        preprocess_fastq_files('data/' + session_id + '/' + renamed_fw_file, 'data/' + session_id + '/' + renamed_rc_file, session_id)
                        
                        flash('Files were successfully uploaded!')
                        return redirect(url_for('web.indexNewSession'))
                    except Exception as e:
                        if os.path.exists('data/'+session_id):
                            rmtree('data/'+session_id)
                        session.clear()
                        print(e) # print error on server side
                        flash('An error occurred while parsing the input files, please make sure the '
                              'files conform the fastq standard')
                        return redirect(url_for('web.indexNewSession'))
                else:
                    flash("Files are already uploaded. Please reset the files first")
                    return redirect(url_for('web.indexNewSession'))
        else:
            flash('Unsupported file types')
            return redirect(url_for('web.indexNewSession'))

    return render_template('indexNewSession.html', form=form, scroll = "import")

# Function that allows to upload a TSV file
@bp.route('/indexLoadSession', methods = ['GET','POST'])
def indexLoadSession():
    form = TSVForm()            # This variable is used to create the import buttons and their functionality
    try:
        didBlast = request.args['didBlast']
    except Exception:
        didBlast = None

    # Code that is used when the user press the submit button:
    if form.validate_on_submit():
        
        # Functions to check if the file type is correct:
        tsvfile = secure_filename(form.tsv_file.data.filename)

        # When the files are of the correct filetype use the following code:
        if allowed_file(tsvfile):
            try:
                session_id = get_sessionID()
                if os.path.exists('data/'+session_id):
                    if session.get("id"):
                        session.pop("id")
                    session_id = get_sessionID()
            finally:
                # Rename the files:
                renamed_tsv_file = 'data/' + session_id + '/' + 'tsv_file.'+tsvfile.rsplit('.', 1)[1].lower()
                
                # Create directory and save the files in it:
                
                if not os.path.exists('data/'+session_id):
                    try:
                        os.makedirs('data/'+session_id)
                        form.tsv_file.data.save(renamed_tsv_file)
                        
                        # Create the FastQDataFrame object with the files and export it as a pickle:
                        didBlast = preprocess_tsv_file(renamed_tsv_file, session_id)
                        
                        flash('Files were successfully uploaded!')
                        return redirect(url_for('web.indexLoadSession', didBlast = didBlast))
                    except Exception as e:
                        if os.path.exists('data/'+session_id):
                            rmtree('data/'+session_id)
                        session.clear()
                        print(e) # print error on server side
                        flash(e)
                        return redirect(url_for('web.indexLoadSession'))
                else:
                    flash("Files are already uploaded. Please reset the files first")
                    return redirect(url_for('web.indexNewSession'))
        else:
            flash('Unsupported file types')
            return redirect(url_for('web.indexLoadSession'))

    return render_template('indexLoadSession.html', form=form, scroll = "import", didBlast = didBlast)

def get_sessionID():
    try:
        session_id = session['id']          # Session_id is the key to know which directory the files need to be saved to
    except KeyError:
        session_id = str(uuid.uuid1())      # Create an unique ID code
        session['id'] = session_id
    return session_id

@bp.route('/blastPage', methods=['GET', 'POST'])
def blastPage():
    #get the desired form class
    form = XMLForm()
    #check if the form has anython in it
    if form.validate_on_submit():
        #get the file name
        xml_file = secure_filename(form.xml_file.data.filename)
        #check if the extension is allowed
        if allowed_file(xml_file):
            #try to get the session else make one
            try:
                session_id = session['id']
            except: 
                session_id = str(uuid.uuid1())
                session['id'] = session_id
            finally:
                #check if the path exist and if an pickel exist else make one dir and an empty pickle
                if not os.path.exists('data/'+session_id):
                    os.makedirs('data/'+session_id)
                if not os.path.exists('data/'+session_id+'/'+'pickle.pkl'):
                    initialize_empty_dataframe(session_id)   
                #add the xml file to the session dir
                renamed_fw_file = 'blastout.xml'
                form.xml_file.data.save('data/' + session_id + '/' + renamed_fw_file)
                flash('Files were successfully uploaded!')
                return redirect(url_for('web.blastPage', upload=True, haveResults = False))

        else:
            flash('Unsupported file types')
            return redirect(url_for('web.blastPage', upload=False, haveResults = False))
    df = FastQDataframe.load_pickle("data/" + session['id'] + "/pickle.pkl")
    return render_template('blastPage.html', form=form, upload=False, haveResults = str("fw_full_tax" in df.get_dataframe().columns))


@bp.route('/results', methods=['GET', 'POST'])
def results():
    return render_template('results.html')

@bp.route('/about')
def about():
    return render_template('about.html')
