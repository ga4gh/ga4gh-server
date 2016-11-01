# GA4GH API Example: gemini database

Here is an example of adopting the GA4GH API to an existing and external project: the Gemini software framework for exploring genome variation.

This guide makes use of these instructions: https://github.com/ga4gh/schemas/blob/master/doc/source/appendix/swagger.rst .

## Prerequisites

You'll need swagger-codegen (you'll just need to be able to execute the jar file really: https://github.com/swagger-api/swagger-codegen), the schemas repository (https://github.com/ga4gh/schemas)

git clone https://github.com/ga4gh/schemas.git
git clone https://github.com/swagger-api/swagger-codegen.git
cd swagger-codegen; mvn clean package; cp ./modules/swagger-codegen-cli/target/swagger-codegen-cli.jar ./codegen.jar; java -jar ./codegen.jar help

## Instructions

1. First, we have to take the schemas (\*.proto) defined in the 'schemas' repository, and produce swagger API definitions (\*.swagger.json) from them, via the _protoc_ utility.

    cd schemas
    mkdir -p target/swagger

    Since the swagger-codegen.jar command line utility only accepts 1 API definition file as input, our services will need to all be defined in one file.  Currently, in the schemas repository, they are defined as many files.  I've manually combined them all in to 1 file called _all_services.proto_.  Let's copy that file to the same location as the others:

    cp path/to/all_services.proto ./src/main/proto/ga4gh

    Now we use the _protoc_ utility:

    protoc -I./src/main/proto --swagger_out=. src/main/proto/ga4gh/all_service.proto

    This should produce a directory 'ga4gh', and inside it, all_service.swagger.json ; this is the swagger API service definition file.


2. Now we use the all_service.swagger.json swagger file from step 1 to produce stub server code, via the _swagger-codegen.jar_ utility.

    To support Python 2 in the output stub server, create a JSON configuration file with {"supportPython2":"True"} in it.  Pass this file to swagger-codegen-cli's 'generate' command, via the -c option; we'll call it temp.json .

    java -jar swagger-codegen.jar generate -i schemas/target/swagger/ga4gh/all_service.swagger.json -l python-flask -o /tmp/server -c temp.json

    The folder "server" should then be a Python flask-based stub server that one needs to modify.  You need to make it so that the server implements the GA4GH interface, and calls the Gemini database in the background.

    Rename all of the modules in server/controllers/* so that they go from like_this_controller.py to LikeThis_controller.py .  The reason this step is needed, is because, for some reason, the generated server code expects the modules to look like that, and I don't know how to make it expect otherwise, so this is a workaround.

    mv server/controllers/read_service_controller.py server/controllers/ReadService_controller.py

    Also, go inside each one of those files and change the paramaters to each function from CamelCase to underscore_case.  The server will fail to start without performing this step.  So for example, if inside one of the files you have a definition like "def my_function(datasetId):", change it to "def my_function(dataset_id):".

3. Now we load an example VCF file into the gemini SQLite database:

    gemini load -v file.vcf my.db

    Start the server with:

    python server/app.py

    *Navigate to http://localhost:8080/ui/

    Edit the files under 'controllers' in your stub server.  Their names may need to be modified from something_like_this_controller.py to SomethingLikeThis_controller.py .  Inside each file though, you may have to change variable names from CamelCase to under_score_case
