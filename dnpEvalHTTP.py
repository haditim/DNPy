from http.server import HTTPServer, BaseHTTPRequestHandler
from optparse import OptionParser
import time
import webview  # pip install pywebview && sudo apt-get install python-gi gir1.2-webkit-3.0
import threading


class RequestHandler(BaseHTTPRequestHandler):

    def do_GET(self):

        request_path = self.path
        testVal = time.clock()
        print("\n----- Request Start ----->\n")
        print("Request path:", request_path)
        print("Request headers:", self.headers)
        print("<----- Request End -----\n")
        if 'eval' in request_path:
            message = """
            <html>
            <meta http-equiv="refresh" content="1" />
              <body>
                <h1>Eval Page</h1>
                """ + str(testVal) + """
              </body>
            </html> 
            """
        else:
            message = """
         <html>
         <body>
           <h1>Form Page</h1>
        
           <form method=POST>
           <ul>
            <li>Phasing routine:  
                <select name="phase" size="1">
                  <option value="first">Phase only first dataset</option>
                  <option value="all">Phase all</option>
                  <option value="none">Phase none</option>
                </select> 
            </li>
            <li>T1 calculation data:  
                <select name="t1Calc" size="1">
                  <option value="PCreal">Phase cycled real channel</option>
                  <option value="PCmagn">Phase cycled magnitude channel</option>
                  <option value="real">Real channel</option>
                  <option value="magn">Magnitude channel</option>
                </select> 
            </li>
            <li>Fourier Transform window:
                <input type="number" name="ftWindow">
            
            </li>
            <li>
            <input name="folder" type="file" webkitdirectory directory multiple/>
            </li>
           </ul>
           <input type='submit' value='Submit Form'>
           </form>
        
         </body>
        </html>
        """
        self.send_response(200)
        self.send_header("Set-Cookie", "foo=bar")
        self.end_headers()

        self.wfile.write(bytes(message, "utf8"))

    def do_POST(self):

        request_path = self.path
        self.returnedData = {}
        print("\n----- Request Start ----->\n")
        print("Request path:", request_path)

        request_headers = self.headers
        content_length = request_headers.get('Content-Length')
        length = int(content_length) if content_length else 0

        print("Content Length:", length)
        print("Request headers:", request_headers)
        formData = self.rfile.read(length).decode("utf-8").split('&')
        for i, val in enumerate(formData):
            print(val)
            self.returnedData.update({val.split('=')[0]: val.split('=')[1]})
        print('returnedData:', self.returnedData)
        print("<----- Request End -----\n")
        message1 = """
        <html>
        <body>
   <h1>Form Page</h1>"""
        message2 = """
 </body>
</html>
"""

        # self.send_response(200)
        # self.end_headers()
        # self.wfile.write(bytes(message1+'POST successful'+message2, "utf8"))
        self.send_response(301)
        self.send_header('Location', 'eval')
        self.end_headers()

    do_PUT = do_POST
    do_DELETE = do_GET


def main():
    port = 8080
    print('Listening on localhost:%s' % port)
    server = HTTPServer(('', port), RequestHandler)
    t = threading.Thread(target=server.serve_forever)
    t.start()
    webview.create_window('Evaluation', 'http://localhost:8080', debug=True)
    server.shutdown()


if __name__ == "__main__":
    parser = OptionParser()
    parser.usage = ("Creates an http-server that will echo out any GET or POST parameters\n"
                    "Run:\n\n"
                    "   reflect")
    (options, args) = parser.parse_args()

    main()
