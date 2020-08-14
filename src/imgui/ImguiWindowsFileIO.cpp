#include "ImguiWindowsFileIO.hpp"

#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <filesystem>

#if defined(_WIN32) || defined(_WIN64)
  #include <windows.h>
  #include <direct.h>
  #include <tchar.h>
  #include <IO.h>
  #define GetCurrentDir _getcwd
#else
  #include <unistd.h>
  #include <dirent.h>
  #include <sys/stat.h>
  #define GetCurrentDir getcwd
#endif

#ifdef _MSC_VER 
//not #if defined(_WIN32) || defined(_WIN64) because we have strncasecmp in mingw
  #define strncasecmp _strnicmp
  #define strcasecmp _stricmp
#endif

// is already defined
//#define IM_ARRAYSIZE(_ARR)  ((int)(sizeof(_ARR)/sizeof(*_ARR)))

#if defined(ICON_FA_CARET_DOWN)
  #define CARET_DOWN ICON_FA_CARET_DOWN
#else
  #define CARET_DOWN "v"
#endif

#undef min

using namespace std;
using namespace ImGui;

#if defined(_WIN32) || defined(_WIN64)
vector<string> getWindowsDrives()
{
    vector<string> result;
    TCHAR g_szDrvMsg[] = _T("A:");

    ULONG uDriveMask = _getdrives();

    if (uDriveMask == 0)
    {}
    else
    {
        while (uDriveMask)
        {
            if (uDriveMask & 1)
                result.push_back( g_szDrvMsg );

            ++g_szDrvMsg[0];
            uDriveMask >>= 1;
        }
    }

    return result;
}
#endif 

vector<string> stringSplit(const string &s, char delim)
{
    vector<string> elems;
    string item;

    // use stdlib to tokenize the string
    stringstream ss(s);
    while (getline(ss, item, delim))
        if(!item.empty())
        {
            elems.push_back(item);
        }

    return elems;
}

bool compareNoCase(const string& first, const string& second)
{
    return strcasecmp( first.c_str(), second.c_str() ) <= 0;
}

string MiniPath::filePath() const
{
    return path + name;
}

void MiniPath::fromString (const string& file_path, char delim)
{
    int last_delim_pos = file_path.find_last_of (delim);
    name = file_path.substr (last_delim_pos+1);
    path = file_path.substr (0, last_delim_pos+1);
}

MiniPath::MiniPath()
{}

MiniPath::MiniPath( const string& some_path )
{
    string s = some_path;
    if (MiniPath::isAbsoluteFilePath (s))
    {
        if (s.find ("/") != string::npos) // linux style
            fromString (s ,'/');
        else if (s.find ("\\") != string::npos) // windows style
            fromString (s ,'\\');
    }
    else
    {
        string current = MiniPath::getCurrentDir();

        if (current.find ("/") != string::npos)
        {
            std::replace( s.begin(), s.end(), '\\', '/');
            fromString (current + "/" + s ,'/');
        }
        else if (current.find ("\\") != string::npos)
        {
            std::replace( s.begin(), s.end(), '/', '\\');
            fromString (current + "\\" + s ,'\\');
        }
        else
            fromNameInCurrentDir( s );
    }
}

string MiniPath::getDelim() const
{
    if (path.find ("/") != string::npos) // linux style
        return "/";
    else if (path.find ("\\") != string::npos) // windows style
        return "\\";
    else
    {
        return MiniPath::getSystemDelim();
    }
}

std::string MiniPath::getSystemDelim()
{
    string current = MiniPath::getCurrentDir();
    if (current.find ("/") != string::npos) // linux style
        return "/";
    else if (current.find ("\\") != string::npos) // windows style
        return "\\";
    else
        return "/";
}

string MiniPath::prefix() const
{
    return name.substr (0, name.find_last_of ('.'));
}

string MiniPath::extension() const
{
    if (name.find (".") == string::npos)
        return "";

    return name.substr (name.find_last_of ('.')+1);
}

void MiniPath::addExtension(const string ext) {
   name += '.'+ext;
}

string MiniPath::getCurrentDir()
{
    char cCurrentPath[FILENAME_MAX];
    if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
        return 0;
    cCurrentPath[sizeof(cCurrentPath) - 1] = '\0';
    return string (cCurrentPath);
}

void MiniPath::fromNameInCurrentDir( const string& file_name )
{
    path = getCurrentDir();
    name = file_name;
}

string MiniPath::getName() const
{
    return name;
}

string MiniPath::getPath() const
{
    return path;
}

vector<string> MiniPath::getPathTokens() const
{
    if( !getDelim().empty() )
        return stringSplit( path + getDelim(), getDelim()[0] );

    vector<string> tmp;
    return tmp;
}

void MiniPath::setName( const string& name )
{
    string tmp;
    size_t last_delim_pos = name.find_last_of ('/');
    if( last_delim_pos != string::npos )
        tmp = name.substr (last_delim_pos+1);
    else
        tmp = name;

    size_t last_delim_pos2 = tmp.find_last_of ('\\');
    if( last_delim_pos2 != string::npos )
        tmp = tmp.substr (last_delim_pos2+1);

    this->name = tmp;
}

bool MiniPath::setPath( const string& absolut_path )
{
    if (isAbsoluteFilePath( absolut_path ))
    {
        path = absolut_path;
        return true;
    }
    return false;
}

bool MiniPath::isAbsoluteFilePath( const string& s )
{
    if (s.size() > 0)
    {
        if (s.at(0) == '/' || s.at(1) == ':')
            return true;
        else
            return false;
    }
    else
        return false;
}

std::list<string> MiniPath::listDirectories( const string& s )
{
#if defined(_WIN32) || defined(_WIN64)
	using namespace filesystem;
    list<string> directories;

    struct _finddata_t c_file;
    intptr_t hFile;

    string filter = s + MiniPath::getSystemDelim() + "*.*";

    if( (hFile = _findfirst( filter.c_str(), &c_file )) == -1L )
    {}
    else
    {
        do
        {
            if( ( c_file.attrib & _A_SUBDIR ) )
                directories.push_back( string( c_file.name ) );
        } while( _findnext( hFile, &c_file ) == 0 );
    }


    return directories;
#else

    list<string> directories;

    string path = s;

    DIR* directory = opendir(path.c_str());
    if(directory != NULL)
    {

        struct dirent* entry;
        while( (entry = readdir(directory)) != NULL )
        {
            if(entry->d_type == DT_DIR)
                directories.push_back(string(entry->d_name));
        }

        closedir(directory);
    }

    directories.sort(compareNoCase);

    return directories;
#endif
}

std::list<string> MiniPath::listFiles( const string& s, string filter )
{
#if defined(_WIN32) || defined(_WIN64)
    using namespace filesystem;    list<string> files;


    struct _finddata_t c_file;
    intptr_t hFile;

    string find_filter = s + MiniPath::getSystemDelim() + filter;

    if( (hFile = _findfirst( find_filter.c_str(), &c_file )) == -1L )
    {}
    else
    {
        do
        {
            if( !( c_file.attrib & _A_SUBDIR ) )
                files.push_back( string( c_file.name ) );
        } while( _findnext( hFile, &c_file ) == 0 );
    }


    return files;
#else

    list<string> directories;

    string path = s;

    DIR* directory = opendir(path.c_str());
    if(directory != NULL)
    {

        struct dirent* entry;
        while( (entry = readdir(directory)) != NULL )
        {
            if(entry->d_type == DT_REG)
            {
                vector<string> splits = stringSplit(string(entry->d_name),'.');
                vector<string> filtersplit = stringSplit(filter, '.');

                if(filter.empty() || filtersplit.back().compare("*") == 0 ||  splits.back().compare(filtersplit.back()) == 0)
                    directories.push_back(string(entry->d_name));
            }
        }

        closedir(directory);
    }

    directories.sort(compareNoCase);

    return directories;
#endif

}

bool MiniPath::pathExists( const std::string& s )
{
#if defined(_WIN32) || defined(_WIN64)
    struct _finddata_t c_file;
    intptr_t hFile;

    if( (hFile = _findfirst( s.c_str(), &c_file )) != -1L )
        return true;

    return false;
#else

    struct stat info;

    if(  stat(s.c_str(), &info) == 0 )
        return true;

    return false;
#endif
}

vector<string> toLimited( 
    const vector<string>& vs, 
    int readable_length = 5 )
{
    vector<string> res( vs.size() );

    size_t count = 0;
    for( const string& s : vs )
    {
        if( readable_length < 5 )
           res[count++] = s.c_str();
        else    
        {
            string tmp2 = s.substr( 0, readable_length / 2 - 1 );
            tmp2 += "...";
            tmp2 += s.substr( s.size() - readable_length / 2 + 1, s.size()-1 );
            res[count++] = tmp2.c_str();
        }
    }

    return res;
}

vector<string> toLimited( 
    const list<string>& vs, 
    int readable_length = 5 )
{
    vector<string> res( vs.size() );

    size_t count = 0;
    for( const string& s : vs )
    {
        if( readable_length < 5 )
           res[count++] = s.c_str();
        else    
        {
            string tmp2 = s.substr( 0, readable_length / 2 - 1 );
            tmp2 += "...";
            tmp2 += s.substr( s.size() - readable_length / 2 + 1, s.size()-1 );
            res[count++] = tmp2.c_str();
        }
    }

    return res;
}

void toCStringVec( 
    vector<const char*>& res, 
    const vector<string>& vs )
{
    res.resize( vs.size() );
    size_t count = 0;

    for( const string& s : vs )
        res[count++] = s.c_str();
}


bool fileIOWindow(
        bool& try_open,
        string& file_path,
        const std::vector<string>& recently_used_files,
        const string& button_text,
        const std::vector<std::string> file_filter,
        bool ensure_file_exists,
        ImVec2 size )
{
    bool close_it = false;

    vector<const char*> extension_cstrings; 
    toCStringVec( extension_cstrings, file_filter );
    
    static string current_folder = "x";
    static char c_current_folder[ 2048 ];

    if( current_folder == "x" )
        current_folder = MiniPath::getCurrentDir();

    static char current_file[ 256 ] = "select_a_file.json";
    static int  file_type_selected = 0;
    static int  file_selected = 0;
    static int  directory_selected = 0;
#if defined(_WIN32) || defined(_WIN64)
    static int  drive_selected = 0;
#endif
    static bool directory_browsing = false;
    static int  recent_selected = 0;

    string sys_delim = MiniPath::getSystemDelim();

    string tmp = current_folder + sys_delim + current_file;
    MiniPath current_mini_path( tmp );
    
    list<string> directories = MiniPath::listDirectories( current_mini_path.getPath() );
    std::vector<const char*> dir_list;
    for( const string& s : directories )
    {
        if( s == "." ) continue;
        dir_list.push_back( s.c_str() );
    }
    
    list<string> local_files = MiniPath::listFiles( current_mini_path.getPath(), string( extension_cstrings[file_type_selected] ) );
    vector<string> local_files_readable = toLimited( local_files, 0 );
    vector<const char*> file_list;

    toCStringVec( file_list, local_files_readable );

    if ( directory_browsing )
        size.y += std::min( size_t( 8 ), std::max( dir_list.size(), file_list.size() ) ) *  GetFontSize();

    SetNextWindowSize( size, ImGuiCond_FirstUseEver);
    // I have no idea why I can't resize this window!?!
    ImGuiWindowFlags window_flags = 0;//!ImGuiWindowFlags_AlwaysAutoResize;
    Begin( "FileIO", NULL, window_flags );

    Text("Directory: "); SameLine();
    PushItemWidth( GetWindowWidth() - 145 );

    strcpy( c_current_folder, current_folder.c_str() );
    if( InputText( " ", c_current_folder, IM_ARRAYSIZE( c_current_folder ) ) )
    {
        string tmp = c_current_folder;
        if( MiniPath::pathExists( tmp ) )
            current_folder = c_current_folder;
    }

    if( !recently_used_files.empty() )
    {
    	SameLine();

        std::vector<string> recent_readable = toLimited( recently_used_files, 27 );
        std::vector<const char*> recent;
        toCStringVec( recent, recent_readable );

        if( Button(" " CARET_DOWN " ") )
            ImGui::OpenPopup("RecentFiles");

        if (ImGui::BeginPopup("RecentFiles"))
        {
            if( ListBox( "", &recent_selected, recent.data(), recent.size() ) )
            {
                current_mini_path.fromString( recently_used_files[recent_selected], sys_delim[0] );
                if( !current_mini_path.getName().empty() )
                    strcpy( current_file, current_mini_path.getName().c_str() );
                current_folder = current_mini_path.getPath();
                CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }
    }
    Text("           "); SameLine();
    vector<string> split_directories = current_mini_path.getPathTokens();
    for( size_t i = 0; i < split_directories.size(); ++i )
    {
        if( Button( split_directories[i].c_str() ) )
        {
#if defined(_WIN32) || defined(_WIN64)
            string chosen_dir;
#else
            string chosen_dir = "/";
#endif

            for( size_t j = 0; j <= i; ++j)
            {
                chosen_dir += split_directories[j];
                if( j != i )
                    chosen_dir += sys_delim;
            }
            current_folder = chosen_dir;
        }
        if( i != split_directories.size()-1 )
            SameLine();
    }

    if( CollapsingHeader("Browse Directories", ImGuiTreeNodeFlags_DefaultOpen) )
    {
        directory_browsing = true;

#if defined(_WIN32) || defined(_WIN64)
        vector<const char*> drive_list;
        toCStringVec( drive_list, getWindowsDrives() );

        // select list item dependent on folder path
        if( drive_list[drive_selected][0] != c_current_folder[0]  )
            for( size_t i = 0; i < drive_list.size(); ++i )
                if( c_current_folder[0] == drive_list[i][0] )
                    drive_selected = i;

        Text( "  " ); SameLine();
        PushItemWidth( 40 );
        if( ListBox( "  ", &drive_selected, drive_list.data(), drive_list.size() ) )
        {
            string new_path = drive_list[drive_selected];
            current_folder = new_path;
        }
#else
        Text( "           " );
#endif
        SameLine();

        PushItemWidth( (GetWindowWidth()-60)*0.4 );
        if( ListBox( " ", &directory_selected, dir_list.data(), dir_list.size() ) )
        {
            string new_path;

            if( string(dir_list[directory_selected]).find( ".." ) != std::string::npos )
            {
                vector<string> split_directories = current_mini_path.getPathTokens();
                current_folder.clear();
#if defined(_WIN32) || defined(_WIN64)
                // noop
#else
                current_folder += "/";
#endif
                for( size_t i = 0; i < split_directories.size() - 1; ++i )
                {
                    current_folder += split_directories[i];
                    if( i < split_directories.size() - 2 )
                        current_folder += sys_delim;
                }
            }
            else
                current_folder = current_folder + sys_delim + dir_list[directory_selected];
        }

        SameLine();
        PushItemWidth( (GetWindowWidth()-60)*0.6 );
        if( ListBox( "", &file_selected, file_list.data(), file_list.size() ) )
        {
            strcpy( current_file, file_list[file_selected] );
        }
    }
    else
        directory_browsing = false;

    Text(" ");

    Text("File Name: "); SameLine();
    InputText( "  ", current_file, IM_ARRAYSIZE( current_file ) );

    Text( "File Type: " ); SameLine(); Text( "%s", extension_cstrings[file_type_selected] ); SameLine();
    if (Button( CARET_DOWN ) )
        ImGui::OpenPopup("FileType");

    if (ImGui::BeginPopup("FileType"))
    {
        if( ListBox( "", &file_type_selected, extension_cstrings.data(), extension_cstrings.size() ) )
        {
            CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }

    ImVec2 pos = GetCursorPos();

    pos.x += GetWindowWidth() - 75 - button_text.size() * ( GetFontSize() * 0.7 );
    pos.y  = GetWindowHeight() - 30;

    SetCursorPos( pos );
    if( Button( button_text.c_str() ) )
    {
        vector<string> ext_filter_v = stringSplit( extension_cstrings[file_type_selected], '.');
        if (current_mini_path.extension() != ext_filter_v[1]) { current_mini_path.addExtension(ext_filter_v[1]); }
        file_path = current_mini_path.filePath();
        std::cout << file_path << std::endl;
        if( ensure_file_exists )
        {
            if( MiniPath::pathExists( file_path ) )
                close_it = true;
        } else {
            close_it = true;
        }

        try_open = true;
    }

    SameLine();
    if( Button( "Cancel" ) ) {
        close_it = true;
        try_open = false;
    }
    End();

    return close_it;
}
